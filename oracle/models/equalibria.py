#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" Equialibria (excitation & ionisation balance) model for stellar spectra. """

from __future__ import absolute_import, print_function

__author__ = "Andy Casey <arc@ast.cam.ac.uk>"

import logging
import numpy as np
from time import time
from scipy import stats, sparse, ndimage, optimize as op
from astropy import (modeling, table, units as u)

from oracle import (atmospheres, specutils, synthesis, utils)
from oracle.models import Model, utils

logger = logging.getLogger("oracle")

import matplotlib.pyplot as plt


class Converged(BaseException):
    pass

def equalibrium_state(transitions, log_eps, metallicity,
    excitation_regression_species=None, ionisation_state_species=None,
    abundance_state_species=None, rew_regression_species=None):
    """
    Calculate the equalibrium state for the atomic transitions and abundances
    provided.
    """

    finite = np.isfinite(log_eps)

    def match_species(species):
        if species is None:
            return np.ones(len(transitions), dtype=bool)
        return np.array([each in species for each in transitions["species"]])

    # Scale the log_eps abundances to the expected relative values.
    expected_log_eps \
        = metallicity + atmospheres.solar_abundance(transitions["species"])

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
    abundance_state = np.nanmedian(log_eps_differences[mask])

    # Calculate the ionisation state.
    # Each set of single and ionised species gives us an indication to what the
    # ionisation state of the system is. When we have many Fe lines but only 1
    # Ti 1 + 2 set, or vice versa, we should average each set.
    species = set(transitions["species"])
    if ionisation_state_species is not None:
        species = species.intersection(abundance_state_species)

    assert len(species) > 1 # Need at least X I and X II

    # For each species and ionisation state, calculate the abundance.
    species_abundances = {}
    for each in sorted(species):
        mask = match_species([each]) * finite
        values = log_eps_differences[mask]

        N = np.isfinite(values).sum()
        species_abundances[each] = [np.nanmedian(values), N]

    # If we have many different elements, create some weighted average.
    elements = set(map(int, species))   
    multiple_ionising_elements = len(elements) > 1 
    if multiple_ionising_elements:
        logger.debug("Calculating relative weight for ionisation state.")

    relative_weights = np.zeros(len(elements))
    element_ionisation_state = np.zeros(len(elements))
    for i, element in enumerate(elements):
        neutral_log_eps, neutral_N = species_abundances.get(
            float(element), (np.nan, 0))
        ionised_log_eps, ionised_N = species_abundances.get(
            float(element) + 0.1, (np.nan, 0))

        # Weight by the least number of lines (e.g., either ionised or neutral).
        relative_weights[i] = min([neutral_N, ionised_N])
        element_ionisation_state[i] = neutral_log_eps - ionised_log_eps

        if multiple_ionising_elements:
            logger.debug("Ionisation state of {0} is {1:.3f} dex ({2} neutral "\
                "{3} ionised lines)".format(element, element_ionisation_state[i],
                neutral_N, ionised_N))

    relative_weights /= relative_weights.sum()
    ionisation_state = (element_ionisation_state * relative_weights).sum() \
        / len(elements)
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



class EqualibriaModel(Model):

    """
    This class performs excitation and ionization balances in order to
    provide a point estimate of the stellar parameters.
    """

    # Default configuration for a EqualibriaModel class
    _default_config = {
        "model": {
            "redshift": True,
            "instrumental_resolution": True,
            "continuum": False,
            "profile_function": "gaussian",
            "atmosphere": {
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

        super(EqualibriaModel, self).__init__(configuration)

        # Initialise the atomic transitions.
        #self.atomic_transitions = self._initialise_atomic_transitions()

        # Initiate an atmosphere interpolator.
        self._atmosphere_interpolator = atmospheres.interpolator(
            **self.config["model"].get("atmosphere", {}))
        return None


    # For pickling and unpickling the class.
    def __getstate__(self):
        allowed_keys = ("config", "atomic_transitions", "_initial_theta")
        state = self.__dict__.copy()
        for key in state.keys():
            if key not in allowed_keys:
                del state[key]
        return state

    def __setstate__(self, state):
        self.__dict__ = state.copy()



    def estimate_stellar_parameters(self, spectra=None, transitions=None,
        initial_theta=None, transition_solver=None, sigma_clip=2.0, clips=1,
        state_tolerance=1e-8, fixed=None, full_output=False, **kwargs):
        """
        Can provide either data=[spectra] or measured atomic transitions
        """

        # if it's data=[spectra], then the initial_theta may need to include vrad, continuum etc,

        # if it's data=[transitions], then we can really start from anywhere, because it is just teff, logg.

        if spectra is None and transitions is None:
            raise ValueError("give at least spectra or transitions")

        if spectra is not None:
            raise NotImplementedError


        assert transitions is not None
        assert fixed is None

        if initial_theta is None:
            # Assume whatever value.
            initial_theta = {
                "effective_temperature": 5750,
                "surface_gravity": 4.5,
                "metallicity": 0.,
                "xi": 1.0
            }

            initial_theta = [5750, 1.0, 4.5, 0.]
            # original:
            initial_theta = [4500, 1.5, 2.0, -1.5]
            initial_theta = [5750, 1.0, 4.5, 0.]
            
        debug = kwargs.pop("debug", False)
        equalibrium_state_kwds = kwargs.pop("equalibrium_state", {})

        global acceptable, sampled_theta, sampled_state_sums

        sampled_theta = []
        sampled_state_sums = []
        acceptable = np.ones(len(transitions), dtype=bool)

        def objective_function(theta, full_output=False):

            global acceptable, sampled_theta, sampled_state_sums

            logger.debug("Stellar parameters: {}".format(theta))

            _exception_response = np.nan * np.ones(len(theta))
            _exception_full_response = (_exception_response,
                np.nan * acceptable.sum(), {})

            #effective_temperature, surface_gravity, metallicity, xi = theta
            effective_temperature, xi, surface_gravity, metallicity =  theta

            if np.any(~np.isfinite(theta)) \
            or not (8000 >= effective_temperature >= 3000) \
            or 0 > xi:
                return _exception_full_response \
                    if full_output else _exception_response

            try:
                photosphere = self._atmosphere_interpolator(
                    effective_temperature, surface_gravity, metallicity)
                atomic_abundances = synthesis.moog.atomic_abundances(
                    transitions[acceptable], photosphere, microturbulence=xi,
                    debug=debug)

            except:
                state, atomic_abundances, info = _exception_full_response
                logger.exception("Exception while calculating abundances at {}"\
                    .format(theta))

            else:
                # Calculate the slopes w.r.t. excitation potential and REW, etc.
                state, info = equalibrium_state(transitions[acceptable],
                    atomic_abundances, metallicity, **equalibrium_state_kwds)

            # Append to the sample.
            total_state = (state**2).sum()
            sampled_theta.append(np.copy(theta))
            sampled_state_sums.append(total_state)

            logger.debug("Equalibrium state: {0} {1:.3e}".format(
                state, total_state))

            if total_state < state_tolerance and not full_output:
                raise Converged(theta, state, total_state)

            if full_output:
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
                    atmospheres.solar_abundance(
                        transitions["species"][acceptable])

                # Remove anything more than |sigma_clip| from the median.
                sigma \
                    = (log_eps_differences - np.nanmedian(log_eps_differences))\
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
        all_abundances = np.nan * np.ones(len(transitions))
        all_abundances[acceptable] = final_abundances

        results_table = table.Table(transitions.copy())
        results_table.add_column(table.Column(
            name="log_eps", data=np.round(all_abundances, 3)))
        results_table.add_column(table.Column(
            name="outlier", data=~acceptable, dtype=bool))

        if full_output:
            sampled_theta = np.array(sampled_theta)
            sampled_state_sums = np.array(sampled_state_sums)

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
                photosphere = self._atmosphere_interpolator(*stellar_parameters)
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


def wavelengths_in_data(wavelengths, data):
    orders = -np.ones(len(wavelengths), dtype=int)
    for i, spectrum in enumerate(data):
        orders[(spectrum.disp[-1] >= wavelengths) \
            * (wavelengths >= spectrum.disp[0])] = i
    return orders


