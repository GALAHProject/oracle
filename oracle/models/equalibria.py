#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" Equialibria (excitation & ionisation balance) model for stellar spectra. """

from __future__ import absolute_import, print_function

__author__ = "Andy Casey <arc@ast.cam.ac.uk>"

__all__ = ["BaseEqualibriumModel", "EqualibriumModel"]

import logging
import numpy as np
from time import time
from scipy import stats, sparse, ndimage, optimize as op
from astropy import (modeling, table, units as u)

from oracle import (photospheres, solvers, specutils, synthesis, utils)
from oracle.transitions import AtomicTransition
from oracle.models import Model, utils

logger = logging.getLogger("oracle")

import matplotlib.pyplot as plt



class BaseEqualibriumModel(Model):

    """
    This class performs excitation and ionization balances in order to
    provide a point estimate of the stellar parameters.
    """

    # Defaults.
    DEFAULT_TRANSITION_LIMITS = {
        "equivalent_width": [20, 120]
    }
    DEFAULT_INITIAL_THETA = {
        "effective_temperature": 5750,
        "surface_gravity": 4.5,
        "metallicity": 0.,
        "xi": 1.0
    }
    DEFAULT_INITIAL_THETA = [5750, 4.5, 0, 1.0]
    DEFAULT_PARAMETER_LIMITS = {
        "effective_temperature": [3000, 8000],
        "surface_gravity": [0, 5],
        "xi": [0, 5],
        "metallicity": [-5, 0.5]
    }

    # Default configuration for a BaseEqualibriumModel class
    _default_config = {
        "model": {
            "redshift": True,
            "instrumental_resolution": True,
            "continuum": False,
            "profile_function": "gaussian",
            "photosphere": {
                "kind": "marcs"
            }
        },
        "settings": {
            "threads": 1
        }
    }

    def __init__(self, configuration=None, **kwargs):
        """
        Initialise an Equialibria Model class.

        :param configuration:
            The path to the configuration file that contains settings for the
            class, a dictionary describing the configuration, or a string
            containing the configuration in a YAML-dump format. Only YAML-style
            formats are accepted.

        :type configuration:
            str or dict
        """

        super(BaseEqualibriumModel, self).__init__(configuration)

        # Initialise the atomic transitions.
        #self.atomic_transitions = self._initialise_atomic_transitions()

        # Initiate an photosphere interpolator.
        self._photosphere_interpolator = photospheres.interpolator(
            **self.config["model"].get("photosphere", {}))

        self._equalibrium_estimate_ = None

        return None


    # For pickling and unpickling the class.
    def __getstate__(self):
        allowed_keys = ("config", "_equalibrium_estimate_", "_initial_theta")
        state = self.__dict__.copy()
        for key in state.keys():
            if key not in allowed_keys:
                del state[key]
        return state

    def __setstate__(self, state):
        self.__dict__ = state.copy()


    def estimate_stellar_parameters(self,
        atomic_transitions, initial_theta=DEFAULT_INITIAL_THETA,
        transition_limits=DEFAULT_TRANSITION_LIMITS, sigma_clip=2.0, clips=1,
        fixed_parameters=None, parameter_limits=DEFAULT_PARAMETER_LIMITS,
        state_tolerance=1e-8, full_output=False, **kwargs):
        """
        Estimate the stellar parameters by performing excitation and ionisation
        equalibrium.

        This method can operate using an atomic transition table with measured
        equivalent widths, or it can measure atomic transitions from spectra.
        """

        # If transitions is given and includes equivalent_widths, then that's
        # all we need.

        # If spectra is given, then transitions is needed (measured or not),
        # as well as a transition_solver. The transition_solver might also have
        # something like background_transitions.


        # atomic_transitions should have measured equivalent widths!
        measured = np.isfinite(atomic_transitions["equivalent_width"]) \
            * (atomic_transitions["equivalent_width"] > 0)
        if not any(measured):
            raise ValueError("atomic transitions table must have measured "\
                "equivalent widths, or spectra and a transition_solver are"\
                " required")

        assert fixed_parameters is None

        initial_theta = initial_theta.copy()
        transition_limits = transition_limits.copy()
        parameter_limits = parameter_limits.copy()

        # Transition limits can be on wavelength, equivalent_width, rew,
        # log_eps,
        available_transition_limits = ("wavelength", "equivalent_width",
            "reduced_equivalent_width", "log_eps")
        for key in [] + transition_limits.keys():
            if key not in available_transition_limits:
                logger.warn("Ignoring transition limit on unrecognised "
                    "parameter '{}'".format(key))
            del transition_limits[key]

        logger.info("Limiting transitions to range: {}".format(transition_limits))
        logger.info("Limiting parameters to range: {}".format(parameter_limits))

        debug = kwargs.pop("debug", False)
        equalibrium_state_kwds = kwargs.pop("equalibrium_state", {})

        global acceptable, sampled_theta, sampled_state_sums

        sampled_theta = []
        sampled_state_sums = []
        acceptable = np.ones(len(atomic_transitions), dtype=bool)

        # Any transition limits that are not dependent on abundances can be
        # applied now, before the objective function.
        for key, (upper, lower) in transition_limits.items():
            is_ok = (upper >= atomic_transitions[key]) \
                * (atomic_transitions[key] >= lower)
            if not np.any(~is_ok):
                logger.info("{0} transitions ignored due to restriction on {1}"\
                    .format((~is_ok).sum(), key))
            acceptable *= is_ok


        def objective_function(theta, full_output=False):

            global acceptable, sampled_theta, sampled_state_sums

            # Total disaster recovery:
            if (len(sampled_state_sums) > 100 \
                and not np.any(np.isfinite(sampled_state_sums[-100:]))):
                raise ValueError("last hundred sampled thetas returned NaNs")

            _exception_response = np.nan * np.ones(len(theta))
            _exception_full_response = (_exception_response,
                np.nan * acceptable.sum(), {})
                
            def invalid_value():
                sampled_theta.append(theta)
                sampled_state_sums.append(np.nan)

                return _exception_full_response \
                    if full_output else _exception_response 

            logger.debug("Stellar parameters: {}".format(theta))
            if np.any(~np.isfinite(theta)):
                return invalid_value()

            names = \
                ("effective_temperature", "xi", "surface_gravity", "metallicity")
            effective_temperature, xi, surface_gravity, metallicity = theta

            for name, _ in zip(names, theta):
                if name in parameter_limits:
                    lower, upper = parameter_limits[name]
                    if not (upper >= _ >= lower):
                        return invalid_value()

            try:
                photosphere = self._photosphere_interpolator(
                    effective_temperature, surface_gravity, metallicity)
                atomic_abundances = synthesis.moog.atomic_abundances(
                    atomic_transitions[acceptable], photosphere,
                    microturbulence=xi, debug=debug)

            except:
                state, atomic_abundances, info = _exception_full_response
                logger.exception("Exception while calculating abundances at {}"\
                    .format(theta))

            else:

                # Apply any transition limits that depend on abundance.
                # Currently this is just log_eps (#TODO this will need updating)
                if "log_eps" in transition_limits:
                    l, u = transition_limits["log_eps"]
                    still_acceptable = \
                        (u >= atomic_abundances) * (atomic_abundances >= l)
                    if np.any(~still_acceptable):
                        logger.info("{0} transitions ignored at these stellar "\
                            "parameters due to restriction on log_eps".format(
                                (~still_acceptable).sum()))
                else:
                    still_acceptable = np.ones(acceptable.sum(), dtype=bool)

                # Calculate the slopes w.r.t. excitation potential and REW, etc.
                state, info = equalibrium_state(
                    atomic_transitions[acceptable][still_acceptable],
                    atomic_abundances[still_acceptable], metallicity,
                    **equalibrium_state_kwds)

            # Append to the sample.
            total_state = (state**2).sum()
            sampled_theta.append(np.copy(theta))
            sampled_state_sums.append(total_state)

            logger.debug("Equalibrium state: {0} {1:.3e}".format(
                state, total_state))

            if total_state < state_tolerance and not full_output:
                raise Converged(theta, state, total_state)

            if full_output:
                # NOTE: atomic_abundances here is only for *acceptable* lines.
                return (state, atomic_abundances, info)
            return state

        # Optimisation keywords (fsolve and fmin)
        op_fsolve_kwds = {
            "fprime": utils.jacobian_original,
            "col_deriv": 1,
            "epsfcn": 0,
            "xtol": 1e-10,
            "maxfev": 100
        }
        op_fsolve_kwds.update(kwargs.pop("op_fsolve_kwargs", {}))
        op_fsolve_kwds["full_output"] = True

        op_fmin_kwds = {
            "xtol": 0.1,
            "ftol": 0.1,
            "maxiter": 500,
            "maxfun": 500,
            "disp": False
        }
        op_fmin_kwds.update(kwargs.pop("op_fmin_kwargs", {}))
        op_fmin_kwds["full_output"] = True

        _exception_response = (np.nan * np.ones(4), None, np.nan * np.ones(4))
        _exception_full_response = (np.nan * np.ones(4), None,
            np.nan * np.ones(4), sampled_theta, sampled_state_sums, {})

        iteration, t_init = 0, time()
        
        while True:

            try:
                x, info_dict, ier, mesg = op.fsolve(
                    objective_function, initial_theta, **op_fsolve_kwds)

            except Converged as e:
                logger.info("Stopped optimisation function because convergence "
                    "has been absolutely achieved.")

                # e.args contains (theta, state, total_state)
                x = np.array(e.args[0])
                ier, mesg = 1, "Optimisation absolutely converged."
                info_dict = dict(zip(
                    ("nfev", "njev", "fvec", "fjac", "r", "qtf"), [np.nan] * 6))

            except ValueError:
                logger.exception("Unrecoverable exception occurred during the "
                    "optimisation:")
                return _exception_full_response \
                    if full_output else _exception_response

            # If the fsolve optimisation fails, we should try again.
            if ier != 1:
                logger.warn(
                    "Initial optimisation failed with the message below. Now "\
                    "attempting a dumber, more robust optimisation. Error: {}"\
                    .format(mesg))

                y = lambda x, **k: (objective_function(x, **k)**2).sum()
                try:
                    x, fopt, n_iter, funcalls, warnflag = op.fmin(y, x, 
                        **op_fmin_kwds)

                except Converged as e:
                    logger.info("Stopped optimisation function because "
                        "convergence has been absolutely achieved.")
                    # e.args contains (theta, state, total_state)
                    x = np.array(e.args[0])
                    fopt = np.array(e.args[2])
                    n_iter, funcalls, warnflag = -1, -1, 0

                except ValueError:
                    logger.exception("Unrecoverable exception occurred during "
                        "the optimisation:")
                    return _exception_full_response \
                        if full_output else _exception_response

                if warnflag != 0:
                    logger.warn([
                        "Maximum number of evaluations made by Nelder-Mead.",
                        "Maximum number of iterations reached by Nelder-Mead."
                        ][warnflag - 1])
                else:
                    logger.info("Nelder-Mead optimisation completed successfully")
            else:
                fopt, n_iter, funcalls, warnflag = np.nan, 0, 0, -1
                logger.info("fsolve optimisation completed successfully")

            # Update the number of iterations.
            iteration += 1

            # Do any outlier clipping.
            if clips >= iteration:

                metallicity = x[3]
                state, log_eps, info = objective_function(x, full_output=True)
                
                # We will fit the regression lines to the scaled-solar difference.
                log_eps_differences = log_eps - metallicity - \
                    photospheres.solar_abundance(
                        atomic_transitions["species"][acceptable])

                # Remove anything more than |sigma_clip| from the mean/median.
                avg = np.nanmean if equalibrium_state_kwds.get(
                    "averaging_behaviour", "mean") else np.nanmedian
                sigma \
                    = (log_eps_differences - avg(log_eps_differences))\
                        / np.nanstd(log_eps_differences)

                outlier = np.abs(sigma) > sigma_clip
                logger.info("On iteration {0}, {1} outlier(s) were removed"\
                    .format(iteration, outlier.sum()))
                acceptable[acceptable] *= ~outlier

                # Update the initial theta with the optimised result.
                initial_theta = x.copy()

                if not np.any(outlier):
                    logger.info("Optimisation complete.")
                    break

            else:
                logger.info("Optimisation complete.")
                break

        optimisation_info = {
            "fsolve": {
                "nfev": info_dict["nfev"],
                "njev": info_dict["njev"],
                "fvec": info_dict["fvec"],
                "fjac": info_dict["fjac"],
                "r":    info_dict["r"],
                "qtf":  info_dict["qtf"],
                "ier":  ier,
                "mesg": mesg,
                "kwds": op_fsolve_kwds
            },
            "fmin": {
                "fopt": fopt,
                "niter": n_iter,
                "funcalls": funcalls,
                "warnflag": warnflag,
                "kwds": op_fmin_kwds
            },
            "clipping_iterations": iteration,
            "time_taken": time() - t_init
        }

        # Note, remember final_abundances will have size of sum(acceptable) !!
        final_state, final_abundances, info = objective_function(x, True)

        # Create an abundance *table*
        all_abundances = np.nan * np.ones(len(atomic_transitions))
        all_abundances[acceptable] = final_abundances

        results_table = table.Table(atomic_transitions.copy())
        results_table.add_column(table.Column(
            name="log_eps", data=np.round(all_abundances, 3)))
        results_table.add_column(table.Column(name="log_eps_Solar",
            data=photospheres.solar_abundance(results_table["species"])))

        # x = (teff, xi, logg, metallicity)
        results_table.add_column(table.Column(
            name="[X/M]", data=np.round(all_abundances - x[3] \
                - photospheres.solar_abundance(atomic_transitions["species"]), 3)))
        results_table.add_column(table.Column(
            name="outlier", data=~acceptable, dtype=bool))

        # Arrayify for later.
        sampled_theta = np.array(sampled_theta)
        sampled_state_sums = np.array(sampled_state_sums)

        # Save the successful equalibrium information for pickling.
        logger.info("Saving successful equalibrium information to the model.")
        self._equalibrium_estimate_ = (x, results_table, final_state,
            sampled_theta, sampled_state_sums, optimisation_info) 

        if full_output:
            
            # TODO: Return profile fits + results from line fitting solver, where
            # applicable
            return (x, results_table, final_state, sampled_theta,
                sampled_state_sums, optimisation_info)

        return (x, results_table, final_state)



    def estimate_stellar_parameters_old(self, data, initial_theta=None,
        full_output=False, **kwargs):
        """
        Return point estimates for the stellar parameters (effective temperature,
        surface gravity, metallicity, microturbulence) using an excitation and
        ionisation balance approach.

        :param data:
            The observed data.

        :type data:
            list of `oracle.specutils.Spectrum1D` objects

        :param initial_theta: [optional]
            Initial estimates of the stellar parameters, radial velocities, and
            continuum parameters (where appropriate).

        :type initial_theta:
            dict

        :returns:
            Point estimates of the effective temperature, surface gravity,
            metallicity, and microturbulence.
        """

        # TODO: this is just for debugging.
        plot = kwargs.pop("plot", False)
        if initial_theta is None:
            initial_theta = self.initial_theta(data)

        # State keywords for later on.
        state_kwds = kwargs.pop("state_kwargs", {"sigma_clip": 5})
        state_kwds["full_output"] = True

        sp_initial_theta = [initial_theta[p] for p in self._stellar_parameters]

        # Were profiles measured, and where were they measured?
        profiles_measured = self.atomic_transitions.meta.get(
            "profiles_given_stellar_parameters", False)
        if not profiles_measured or sp_initial_theta != profiles_measured:
            # Measure them again.
            #fitted_profiles = self.fit_atomic_transitions2(data, **initial_theta)
            fitted_profiles = self.fit_all_atomic_transitions(data, **initial_theta)

        # Create a copy of the transitions for our little 'objective adventure'.
        global transitions
        transitions = self.atomic_transitions.copy()
        for_equalibria = np.isfinite(transitions["equivalent_width"]) \
            * (transitions["equivalent_width"] > 0) \
            * ((transitions["species"] == 26.0) + (transitions["species"] == 26.1))
        # TODO permit elements other than Fe.
        transitions = transitions[for_equalibria]

        # Measure the initial state and record them.
        initial_state, info = utils.equalibrium_state(transitions,
            metallicity=sp_initial_theta[2], **state_kwds)

        transitions = transitions[info["~outliers"]]

        def objective_function(theta):
            """ Minimise the simultaeous equalibrium constraints. """

            global transitions
            stellar_parameters, microturbulence = theta[:3], theta[3]
            try:
                photosphere = self._photosphere_interpolator(*stellar_parameters)
                abundances = \
                    synthesis.moog.atomic_abundances(transitions,
                        photosphere, microturbulence=microturbulence)

            except:
                logger.exception("Exception while calculating abundances at {}"\
                    .format(theta))
                return [np.inf] * len(theta)

            transitions["abundance"] = abundances
            state, info = utils.equalibrium_state(transitions,
                    metallicity=stellar_parameters[2], **state_kwds)
            # Remove outliers from future iterations
            transitions = transitions[info["~outliers"]]

            # TODO: Save transitions table for future introspection.

            # Remove outliers for future fits.
            logger.debug("State at {0}: {1} --> {2:.2e}".format(theta, state,
                (state**2).sum()))
            return state

        # Optimisation
        op_kwds = {
            "fprime": utils.jacobian_original,
            "col_deriv": 1,
            "epsfcn": 0,
            "xtol": 1e-10,
            "maxfev": 100
        }
        op_kwds.update(kwargs.pop("op_kwargs", {}))
        op_kwds["full_output"] = True

        result = op.fsolve(objective_function, sp_initial_theta, **op_kwds)

        # Show the new lines
        import matplotlib.pyplot as plt
        fig = plt.gcf()
        ax = fig.axes

        state, info = utils.equalibrium_state(
            transitions, metallicity=result[0][2], full_output=True)
        coefficients = info["coefficients"]
        ok = np.isfinite(transitions["abundance"])
        x = transitions["excitation_potential"][ok]
        ax[0].plot(x, np.polyval(coefficients[0], x), c='r')
        x = np.log(transitions["equivalent_width"] \
            / transitions["wavelength"])[ok]
        ax[1].plot(x, np.polyval(coefficients[1], x), c='r')



        stellar_parameters = dict(zip(self._stellar_parameters, result[0]))
        raise a

        if full_output:
            return (stellar_parameters, result)
        return stellar_parameters



def associate_transitions(spectra, transitions):

    spectrum_indices = []
    for i, transition in enumerate(transitions):
        for j, spectrum in enumerate(spectra):
            if transition in spectrum:
                spectrum_indices.append(j)
                break

        else:
            spectrum_indices.append(np.nan)

    return np.array(spectrum_indices)


class EqualibriumModel(BaseEqualibriumModel):

    def estimate_stellar_parameters(self, spectra, atomic_transitions,
        transitions_solver=None, initial_theta=BaseEqualibriumModel.DEFAULT_INITIAL_THETA, all_transitions=None,
        **kwargs):

        # Check the spectra.
        if not isinstance(spectra, (list, tuple)):
            spectra = [spectra]
        for each in spectra:
            if not isinstance(each, specutils.Spectrum1D):
                raise TypeError("'{}' is not a Spectrum1D object".format(each))

        # Check that all atomic transitions.
        if not isinstance(atomic_transitions, (list, tuple)):
            atomic_transitions = [atomic_transitions]
        for each in atomic_transitions:
            if not isinstance(each, AtomicTransition):
                raise TypeError("'{}' is not an AtomicTransition".format(each))

        # Check the transitions_solver.
        if transitions_solver is None:
            transitions_solver = solvers.SynthesisFitter

        if not isinstance(transitions_solver, solvers.BaseFitter):
            raise TypeError("transitions solver is expected to be a sub-class "\
                "of oracle.solvers.BaseFitter")

        synth_reqd = isinstance(transitions_solver, solvers.BaseSynthesisFitter)

        # Associate spectra to each atomic_transition. We will do the fitting
        # for all lines in each spectrum with a separate class. That way if we
        # are fitting all lines simultaneously, the continuum, etc is handled
        # properly.

        spectrum_indices = associate_transitions(spectra, atomic_transitions)


        # We have:
        # - photospheres
        # - radiative transfer code
        # - atomic transitions (line-specific masks, fitting regions, local continuum degree, RV tolerances)
        # - solver (global masks, global/local continuum/resolution behaviour)
        #   the .fit_all() function in the solver should take all atomic transitions for an observed spectrum, and any other relevant things (e.g., a synthesiser factory)
        #   

        # TODO: The profile and synthesis fitters need updating to be able to
        #       take an AtomicTransition, which will have its own settings about
        #       continuum_degree, masks, etc. rather than being instantiated in
        #       the Solver class itself. The Solver just manages behaviour, and
        #       the fit_all() function manages global behaviour for a given chan

        


        debug = kwargs.pop("debug", False)
        photosphere = None
        theta = [] + initial_theta
        #[initial_theta[k] for k in \
        #    ("effective_temperature", "surface_gravity", "metallicity", "xi")]


        while True:

            stellar_parameters, xi = theta[:3], theta[3]

            fitted_transitions = []
            # Fit all transitions on a per-channel basis.
            for i, spectrum in enumerate(spectra):

                channel_transitions = [transition for j, transition in \
                    zip(spectrum_indices, atomic_transitions) if j == i]


                if synth_reqd:
                    # Create a synthesiser factory.
                    photosphere = self._photosphere_interpolator(*stellar_parameters)
                    synthesiser_factory = lambda atomic_number, wavelength_range: \
                        lambda abundance: oracle.synthesis.moog.synthesise(
                            all_transitions, photosphere, wavelength_range,
                            microturbulence=xi, photospheric_abundances=[
                            atomic_number, abundance], debug=debug)

                    args = (synthesiser_factory, )

                else:
                    raise NotImplementedError

                fitted_transitions.extend(transitions_solver.fit(spectrum,
                    channel_transitions, *args, full_output=True))

            # Use the fitted transition abundances to calculate the state.



        # When we have measured the atomic transitions, use the BaseEqualibriumModel
        # to yield the stellar parameters.

        return None




class Converged(BaseException):
    """
    An exception class that is used to alert optimisation algorithms that the
    objective function has already converged on an optimal point.
    """
    pass


def equalibrium_state(transitions, log_eps, metallicity,
    excitation_regression_species=None, ionisation_state_species=None,
    abundance_state_species=None, rew_regression_species=None,
    averaging_behaviour="mean", relative_weighting_behaviour="total"):
    """
    Calculate the equalibrium state for the atomic transitions and abundances
    provided.
    """

    if relative_weighting_behaviour not in ("num_neutral", "num_ionised",
        "total", "minimum_species"):
        raise ValueError("relative weighting behaviour not known")

    if averaging_behaviour not in ("mean", "median"):
        raise ValueError("averaging behaviour not understood")
    else:
        logger.debug("Using {0} averaging.".format(averaging_behaviour))

    avg = np.nanmean if averaging_behaviour == "mean" else np.nanmedian
    finite = np.isfinite(log_eps)

    def match_species(species):
        if species is None:
            return np.ones(len(transitions), dtype=bool)
        return np.array([each in species for each in transitions["species"]])

    # Scale the log_eps abundances to the expected relative values.
    expected_log_eps \
        = metallicity + photospheres.solar_abundance(transitions["species"])

    # We will fit the regression lines to the scaled-solar difference.
    log_eps_differences = log_eps - expected_log_eps

    # Which lines will be used to calculate the slope with excitation potential
    # and abundance?
    mask = match_species(excitation_regression_species) * finite

    # Calculate the slope with excitation potential and abundance.
    exc_slope, exc_offset, exc_r_value, exc_p_value, exc_stderr \
        = stats.linregress(x=transitions["excitation_potential"][mask],
            y=log_eps_differences[mask])

    # Calculate the slope with reduced equivalent width and line abundance.
    mask = match_species(rew_regression_species) * finite
    rew = np.log(transitions["equivalent_width"] / transitions["wavelength"])
    rew_slope, rew_offset, rew_r_value, rew_p_value, rew_stderr \
        = stats.linregress(x=rew[mask], y=log_eps_differences[mask])

    # Calculate the abundance state.
    mask = match_species(abundance_state_species) * finite
    abundance_state = avg(log_eps_differences[mask])

    # Calculate the ionisation state.
    # Each set of single and ionised species gives us an indication to what the
    # ionisation state of the system is. When we have many Fe lines but only 1
    # Ti 1 + 2 set, or vice versa, we should average each set.
    species = set(transitions["species"])
    if ionisation_state_species is not None:
        species = species.intersection(abundance_state_species)

    if 2 > len(species):
        logger.warn("Only one species found ({}); ionisation state is unknown"\
            .format(list(species)[0]))

    #assert len(species) > 1 # Need at least X I and X II

    # For each species and ionisation state, calculate the abundance.
    species_abundances = {}
    for each in sorted(species):
        mask = match_species([each]) * finite
        values = log_eps_differences[mask]

        N = np.isfinite(values).sum()
        species_abundances[each] = [avg(values), N]

    # If we have many different elements, create some weighted average.
    elements = set(map(int, species))   
    multiple_ionising_elements = len(elements) > 1 
    if multiple_ionising_elements:
        logger.debug("Calculating relative weight for ionisation state using {}"
            " number of lines.".format(relative_weighting_behaviour))

    relative_weights = np.zeros(len(elements))
    element_ionisation_state = np.zeros(len(elements))
    for i, element in enumerate(elements):
        neutral_log_eps, neutral_N = species_abundances.get(
            float(element), (np.nan, 0))
        ionised_log_eps, ionised_N = species_abundances.get(
            float(element) + 0.1, (np.nan, 0))

        element_ionisation_state[i] = neutral_log_eps - ionised_log_eps

        if multiple_ionising_elements:
            logger.debug("Ionisation state of {0} is {1:.3f} dex ({2} neutral,"\
                " {3} ionised lines)".format(
                    element, element_ionisation_state[i], neutral_N, ionised_N))

        if neutral_N == 0 or ionised_N == 0:
            # Relative weights default to zero.
            logger.info("Skipping element {} because of missing lines".format(
                element))
            continue

        if relative_weighting_behaviour == "minimum_species":
            # Weight by the least number of lines present in either species.
            relative_weights[i] = min([neutral_N, ionised_N])

        elif relative_weighting_behaviour == "num_ionised":
            # Weight by the number of ionised lines present.
            relative_weights[i] = ionised_N

        elif relative_weighting_behaviour == "num_neutral":
            # Weight by the number of neutral lines present.
            relative_weights[i] = neutral_N

        elif relative_weighting_behaviour == "total":
            # Weight by the total number of lines present.
            relative_weights[i] = neutral_N + ionised_N

    ionisation_state = np.nansum(element_ionisation_state * relative_weights) \
        / relative_weights.sum()
    if multiple_ionising_elements:
        logger.debug("Final ionisation state is {0:.3f} dex (neutral-ionised)"\
            .format(ionisation_state))

    # Calculate the final state.
    state = np.array([
        exc_slope,
        ionisation_state,
        abundance_state,
        rew_slope
    ])

    state = np.array([
        exc_slope,
        rew_slope,
        0.1 * ionisation_state,
        0.1 * abundance_state,
    ])


    # TODO: Include an info dictionary.
    return (state, None)




if __name__ == "__main__":

    import oracle
    data = [
        oracle.specutils.Spectrum1D.load_GALAH("/Users/arc/research/galah/data/iDR1/data/benchmark/18Sco_1.fits", normalised=True, rest=True),
        oracle.specutils.Spectrum1D.load_GALAH("/Users/arc/research/galah/data/iDR1/data/benchmark/18Sco_2.fits", normalised=True, rest=True),
        oracle.specutils.Spectrum1D.load_GALAH("/Users/arc/research/galah/data/iDR1/data/benchmark/18Sco_3.fits", normalised=True, rest=True),
        oracle.specutils.Spectrum1D.load_GALAH("/Users/arc/research/galah/data/iDR1/data/benchmark/18Sco_4.fits", normalised=True, rest=True)
    ]

    from astropy.table import Table
    table_data = Table.read("/Users/arc/codes/oracle/oracle/solvers/galah_lines_150619.csv")
    transitions = []
    for row in table_data:
        transitions.append(AtomicTransition(wavelength=row["wavelength"],
            species=row["species"], e_low=row["excitation_potential"],
            log_gf=row["loggf"]))
        transitions[-1].v_rad_tolerance = 5

    model = EqualibriumModel()
    model.estimate_stellar_parameters(data, transitions,
        transitions_solver=oracle.solvers.SynthesisFitter(global_continuum=False, global_resolution=False),
        all_transitions=table_data)





