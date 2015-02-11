#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" Equialibria (excitation & ionisation balance) model for stellar spectra. """

from __future__ import absolute_import, print_function

__author__ = "Andy Casey <arc@ast.cam.ac.uk>"

import logging
import numpy as np
from scipy import stats, ndimage, optimize as op
from astropy import (modeling, table, units as u)

from oracle import atmospheres, specutils, synthesis, utils
from oracle.models.model import Model
from oracle.models import transitions
from oracle.models.jacobian \
    import approximate as stellar_parameter_jacobian_approximation

logger = logging.getLogger("oracle")

import matplotlib.pyplot as plt



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
            "profile_function": "gaussian"
        },
        "settings": {
            "threads": 1
        }
    }

    def __init__(self, configuration, **kwargs):
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
        self.atomic_transitions = self._initialise_atomic_transitions()

        # Initiate an atmosphere interpolator.
        atmosphere_kwds = self.config["model"].get("atmosphere_kwargs", {})
        self._interpolator = atmospheres.interpolator(**atmosphere_kwds)
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


    def _initialise_atomic_transitions(self):
        """
        Load the atomic transition information from the configuration.
        """

        # Load atomic transitions.
        filename = self.config["model"].get("atomic_transitions_filename", None)
        if filename is not None:
            format = self.config["model"].get("atomic_transitions_format", None)
            try:
                atomic_transitions = \
                    table.Table.read(filename, format=format)
            except:
                raise ValueError("failed to load atomic transitions filename "
                    "from {} -- try specifying `atomic_transitions_format` in "
                    "the model configuration".format(filename))

            # Delete the reference to the filename in the config, because if we
            # don't then we will run into pickling problems.
            del self.config["model"]["atomic_transitions_filename"]

        else:
            # Must be giving them individually.
            transitions = self.config["model"].get("atomic_transitions", None)
            if transitions is None:
                raise ValueError("no atomic transitions found")
            atomic_transitions = table.Table(data=transitions)

            # Delete the reference to the transitions in the config, because if
            # we don't then we will run into pickling problems later.
            del self.config["model"]["atomic_transitions"]

        # Verify that we have the minimum columns required.
        for column in ("wavelength", "species", "excitation_potential", "loggf"):
            if column not in atomic_transitions.dtype.names:
                raise ValueError("missing requried column '{}' in the atomic "
                    "transitions table".format(column))

        if "clean" not in atomic_transitions.dtype.names:
            logger.warn("No 'clean' column found in atomic transitions. All the"
                "transitions are assumed to have no blending by nearby lines.")

        # Create a default table with the additional columns we require.
        N = len(atomic_transitions)
        columns = [
            table.Column(np.ones(N) * np.nan,
                name="equivalent_width", unit="milliAngstrom"),
            table.Column(np.ones(N, dtype=bool),
                name="clean", dtype=bool),
            table.Column(np.ones(N) * np.nan,
                name="abundance"),
            table.Column([""] * N,
                name="custom_mask", dtype="|S32")
            ]
        for column in columns:
            if column.name not in atomic_transitions.columns:
                atomic_transitions.add_column(column)
            else:
                index = atomic_transitions.colnames.index(column.name)
                old_dtype = atomic_transitions.columns[index].dtype
                if column.dtype != old_dtype:
                    new_column = table.Column(
                        data=atomic_transitions[column.name],
                        name=column.name, dtype=column.dtype)

                    atomic_transitions.remove_column(column.name)
                    atomic_transitions.add_column(new_column)

        # Add some units.
        atomic_transitions["wavelength"].unit = u.Angstrom
        atomic_transitions["excitation_potential"].unit = u.eV
        atomic_transitions["equivalent_width"].unit = "milliAngstrom"

        # Check for custom masks.
        for mask in set(atomic_transitions["custom_mask"]).difference({""}):
            if mask not in self.config["model"].get("custom_mask", {}):
                raise ValueError("cannot find custom mask '{}' in the model "
                    "configuration".format(mask))

        atomic_transitions.sort(["wavelength", "species"])
        return atomic_transitions



    def fit_atomic_transitions(self, data, effective_temperature=None,
        surface_gravity=None, metallicity=None, microturbulence=None,
        wavelength_region=2.5, outlier_modeling=True, max_outlier_profiles=5,
        **kwargs):
        """
        Fit absorption profiles to the atomic transitions in this model and
        account for the nearby blends.

        :param data:
            The observed spectra.

        :type data:
            list of :class:`oracle.specutils.Spectrum1D`

        :param effective_temperature: [sometimes optional]
            The effective temperature for the photosphere.

            When there are nearby blending transitions, these lines need to be 
            synthesised in order to accurately measure the actual transition we
            care about. Thus, if any blending (non-clean) transitions are within
            `wavelength_region` of a clean line then this parameter is required.

        :type effective_temperature:
            float

        :param surface_gravity: [sometimes optional]
            The surface gravity for the photosphere.

            When there are nearby blending transitions, these lines need to be 
            synthesised in order to accurately measure the actual transition we
            care about. Thus, if any blending (non-clean) transitions are within
            `wavelength_region` of a clean line then this parameter is required.

        :type surface_gravity:
            float

        :param metallicity: [sometimes optional]
            The scaled-solar metallicity for the photosphere.

            When there are nearby blending transitions, these lines need to be 
            synthesised in order to accurately measure the actual transition we
            care about. Thus, if any blending (non-clean) transitions are within
            `wavelength_region` of a clean line then this parameter is required.

        :type metallicity:
            float

        :param microturbulence: [sometimes optional]
            The microturbulence for the photosphere. This is not required for
            <3D> models.

            When there are nearby blending transitions, these lines need to be 
            synthesised in order to accurately measure the actual transition we
            care about. Thus, if any blending (non-clean) transitions are within
            `wavelength_region` of a clean line then this parameter is required.

        :type microturbulence:
            float

        :param wavelength_region: [optional]
            The +/- region (in Angstroms) around each atomic transition to
            consider.

        :type wavelength_region:
            float

        :param outlier_modeling: [optional]
            Detect inaccurate fits due to blending lines and account for them.
            If enabled, then a poor fit is detected when the centroid of the 
            atomic transition is not within 5 per cent of the data, or if the
            profile FWHM has hit an upper boundary. When this happens the code
            will detect regions that are most discrepant, and attempt to fit
            them with additional profiles (using the same profile FWHM).

        :type outlier_modeling:
            bool

        :param max_outlier_profiles: [optional]
            The maximum number of outlier profiles to add for each atomic
            transition. This keyword argument is ignored if `outlier_modeling`
            is set to `False`.
        """

        oversampling_rate = int(kwargs.pop("oversampling_rate", 4))
        if 1 > oversampling_rate:
            raise ValueError("oversampling rate must be a positive integer")

        # TODO this should really be changed to a limit based on radial velocity.
        wavelength_tolerance = abs(kwargs.pop("wavelength_tolerance", 0))

        # Allow the user to specify bounds on the data.
        # (And here we will set some sensible ones.)
        common_bounds = kwargs.pop("bounds", {
            "stddev": (0, 0.3),
            "amplitude": (0, 1)
        })
        if "wavelength" in common_bounds or "mean" in common_bounds:
            raise ValueError("apply bounds on the profile location through the "
                "wavelength_tolerance keyword argument")

        # Create a handy function to update the compound model properties.
        def _update_compound_model(compound_model, wavelength):

            # Deal with the line we care about first.
            single_model = hasattr(compound_model, "mean")
            transition_mean_key = ["mean_0", "mean"][single_model]
            if wavelength_tolerance > 0:
                compound_model.bounds[transition_mean_key] = (
                    wavelength - wavelength_tolerance,
                    wavelength + wavelength_tolerance
                )
            else:
                compound_model.fixed[transition_mean_key] = True

            for param_name in compound_model.param_names:
                if param_name.startswith("mean_") \
                and param_name != transition_mean_key:
                    # It's an outlier line. Fix the wavelength.
                    compound_model.fixed[param_name] = True

                # Apply common bounds.
                try:
                    prefix, num = param_name.split("_")

                except ValueError:
                    prefix, num = param_name, "0"

                if prefix in common_bounds.keys():
                    compound_model.bounds[param_name] = common_bounds[prefix]

                # Tie the stddevs to the original absorption profile.
                # The stddevs between the absorption transition we care
                # about and the outlier transition are related by:
                # R = lambda_1/delta_lambda_1 = lambda_2/delta_labmda_2
                if prefix == "stddev" and num != "0":
                    compound_model.tied[param_name] = lambda _: _.stddev_0 * \
                        getattr(compound_model, "mean_{}".format(num))/_.mean_0
            return True


        photosphere = None
        # Interpolate a photosphere if we have the information to do so.
        if None not in (effective_temperature, surface_gravity, metallicity):
            photosphere = self._interpolator(
                effective_temperature, surface_gravity, metallicity)

        # Identify clean transitions that are actually in the data.
        # If data_indices is -1 it means the line was not found in the data
        data_indices = \
            wavelengths_in_data(self.atomic_transitions["wavelength"], data)
        measurable = self.atomic_transitions["clean"] * (data_indices > -1)

        # Create additional columns in the atomic_transitions table if needed.
        if "profile_amplitude" not in self.atomic_transitions.dtype.names:
            self.atomic_transitions.add_column(table.Column(
                name="profile_amplitude", 
                data=[np.nan] * len(self.atomic_transitions)))
        if "profile_stddev" not in self.atomic_transitions.dtype.names:
            self.atomic_transitions.add_column(table.Column(
                name="profile_stddev",
                data=[np.nan] * len(self.atomic_transitions)))

        # Create the subset containing the measurable lines.
        transition_indices = np.where(measurable)[0]
        transitions = self.atomic_transitions[measurable]
        data_indices = data_indices[measurable]

        fitted_profiles = []
        for transition_index, transition, data_index \
        in zip(transition_indices, transitions, data_indices):
        
            wavelength = transition["wavelength"]

            # Look for nearby transitions within the wavelength region and
            # ignore this line.
            blending = wavelength_region >= \
                np.abs(self.atomic_transitions["wavelength"] - wavelength) 
            blending[transition_index] = False

            # Slice the data +/- some region.
            spectrum = data[data_index]
            disp_indices = spectrum.disp.searchsorted([
                wavelength - wavelength_region,
                wavelength + wavelength_region
                ]) + [0, 1]
            x = spectrum.disp.__getslice__(*disp_indices)
            y = spectrum.flux.__getslice__(*disp_indices)
            
            # Calculate an initial stddev value based on the x spacing.
            initial_stddev = 2 * 5 * np.diff(x).mean()
            
            # Apply any custom mask.
            # TODO This should just remove the masked pixels from x and y,
            #      because setting them to NaN will break the fitter.

            
            # TODO the amplitude initial guess will have to be udpated in the 
            #      presence of continuum.

            synthesised_spectra = {}
            _ = x.searchsorted(wavelength)
            initial_amplitude = 1.0 - y[_]

            # Any continuum?
            continuum_order = self._continuum_order(data_index)
            if continuum_order > -1:
                # TODO I specify order and astropy uses degree. Switch to degree!
                #profile_init *= modeling.models.Polynomial1D(continuum_order + 1)

                # Set initial estimates of continuum.
                # TODO

                # This will probably fuck up the parameter names.
                raise NotImplementedError

            # Any synthesis?
            if np.any(blending):
                if photosphere is None:
                    raise ValueError("transition at {0:.3f} has blending lines "
                        "within {1:.0f} (so a synthesis approach is needed) but"
                        " not all stellar parameters were given".format(
                            wavelength, wavelength_region))

                # Synthesise a spectrum (with oversampling)
                synth_pixel_size = np.diff(x).mean()/oversampling_rate
                synth_disp, synth_flux = synthesis.moog.synthesise(
                    self.atomic_transitions[blending], photosphere,
                    microturbulence=microturbulence,
                    wavelength_region=[x.min(), x.max()],
                    wavelength_step=synth_pixel_size)

                # Create a custom class that uses the synthesised spectrum.
                class GaussianAbsorption1D(modeling.Fittable1DModel):

                    amplitude = modeling.Parameter(default=1)
                    mean = modeling.Parameter(default=0)
                    stddev = modeling.Parameter(default=1)
                    # TODO see issue 18 on GitHub

                    @staticmethod
                    def evaluate(x, amplitude, mean, stddev):
                        convolved = ndimage.gaussian_filter1d(
                            synth_flux, stddev/synth_pixel_size)
                        sampled = np.interp(x, synth_disp, convolved, 1, 1)
                        return sampled * (1.0 - \
                            modeling.models.Gaussian1D.evaluate(x, amplitude,
                                mean, stddev))

                profile = GaussianAbsorption1D(
                    mean=wavelength, amplitude=initial_amplitude,
                    stddev=initial_stddev)

                # Save the initial profile + synthesised spectra
                # TODO multiply by some continuum
                synth_only_model = GaussianAbsorption1D(
                    mean=wavelength, amplitude=0, stddev=initial_stddev)
                synthesised_spectra["initial_synthesis"] = synth_only_model(x)
                profile_only_model = modeling.models.GaussianAbsorption1D(
                    mean=wavelength, amplitude=initial_amplitude,
                    stddev=initial_stddev)
                synthesised_spectra["initial_profile"] = profile_only_model(x)

            else:
                profile = modeling.models.GaussianAbsorption1D(
                    mean=wavelength, amplitude=initial_amplitude,
                    stddev=initial_stddev)

                # TODO multiply by some continuum
                synthesised_spectra["initial_synthesis"] = np.ones(len(x))
                synthesised_spectra["initial_profile"] = profile(x)

            # Update the bounds, fixed, and tied properties.
            _update_compound_model(profile, wavelength)

            j, fitter = 1, modeling.fitting.LevMarLSQFitter()
            synthesised_spectra["initial_composite"] = profile(x)

            while True:

                # Fit the profile.
                fitted = fitter(profile, x, y)
                
                # Break here if we have no more outlier modeling to do.
                if not outlier_modeling or j > max_outlier_profiles: break

                # Limitingly-high stddev values are good indicators of nearby
                # lines that have not been accounted for.

                # Having a large % difference between the profile and the data
                # at the transition point is another good indicator of nearby
                # lines that have not been accounted for.
                if (j == 1 and fitted.stddev == profile.bounds["stddev"][1])   \
                or (j > 1 and fitted.stddev_0 == profile.bounds["stddev_0"][1])\
                or not (1.05 > y[_]/fitted(x[_]) > 0.95): #absorption is 5% off

                    # Add a(nother) outlier absorption profile to this model at
                    # the location where a blending line is most likely to be.

                    # But ignore locations near existing lines so we don't just
                    # pile up absorption profiles in the same place.
                    existing_means = [getattr(fitted, key) \
                        for key in fitted.param_names if key[:4] == "mean"]

                    difference = abs(fitted(x) - y)
                    for mean in existing_means:
                        __ = np.clip(x.searchsorted([
                            mean - 3 * initial_stddev,
                            mean + 3 * initial_stddev
                        ]) + [0, 1], 0, len(difference))
                        difference.__setslice__(__[0], __[1], 0)

                    most_discrepant = difference.argmax()

                    # TODO continuum will fuck this up too.
                    profile *= modeling.models.GaussianAbsorption1D(
                        mean=x[most_discrepant], stddev=initial_stddev,
                        amplitude=1.0 - y[most_discrepant])

                    # Update the bounds, fixed, and tied properties.
                    _update_compound_model(profile, wavelength)
                    j += 1

                else:
                    # No outlier treatment required, apparently.
                    break

            # Create final copies of things
            synthesised_spectra["fitted_composite"] = fitted(x)

            # Save the final model fit, and the model parameters.
            parameters = dict(zip(fitted.param_names, fitted.parameters))
            fitted_profiles.append((x, y, synthesised_spectra, parameters))

            # Update the equivalent width
            amplitude = \
                parameters.get("amplitude", parameters.get("amplitude_0", None))
            stddev = parameters.get("stddev", parameters.get("stddev_0", None))

            # Integral of Gaussian = amplitude * sigma * sqrt(2 * pi)
            #              (in mA) *= 1000
            self.atomic_transitions["profile_stddev"][transition_index] = stddev
            self.atomic_transitions["profile_amplitude"][transition_index] = \
                amplitude
            self.atomic_transitions["equivalent_width"][transition_index] = \
                1000 * np.sqrt(2 * np.pi) * amplitude * stddev

        sensible = np.isfinite(self.atomic_transitions["equivalent_width"]) \
            * (self.atomic_transitions["equivalent_width"] > 0)
        # Eliminate all other abundances, then re-calculate them.
        self.atomic_transitions["abundance"] = np.nan
        self.atomic_transitions["abundance"][sensible] = \
            synthesis.moog.atomic_abundances(self.atomic_transitions[sensible],
                photosphere, microturbulence=microturbulence)

        # Supply some metadata to the atomic_transitions table
        self.atomic_transitions.meta["profiles_given_stellar_parameters"] = \
            [effective_temperature, surface_gravity, metallicity, microturbulence]
        self.atomic_transitions.meta["abundances_given_stellar_parameters"] = \
            [effective_temperature, surface_gravity, metallicity, microturbulence]

        return fitted_profiles

        # At this point should we consider re-fitting lines that are deviant
        # from the wavelength vs stddev plot
        for i, (x, y, synthesised_spectra, parameters) in enumerate(fitted_profiles):
        
            fig, ax = plt.subplots()
            ax.plot(x,y,c='k')
            ax.plot(x, synthesised_spectra["initial_profile"], "r", lw=1.5, label="Initial profile")
            ax.plot(x, synthesised_spectra["initial_synthesis"], "b", lw=1.5, label="Initial synthesis")
            ax.plot(x, synthesised_spectra["initial_composite"], "m", label="Initial composite")
            ax.plot(x, synthesised_spectra["fitted_composite"], "g", label="Fitted composite")

            ax.legend()
            for p, v in parameters.items():
                if p[:4] == "mean" and p not in ("mean", "mean_0"):
                    ax.axvline(v, c="r")

        # Only update those with good quality constraints.
        fig, ax = plt.subplots(3)

        # Add profile_* columns to the atomic_transitions table:
        # profile_stddev, profile_amplitude

        ax[0].scatter(self.atomic_transitions["excitation_potential"],
            self.atomic_transitions["abundance"], facecolor='k')
        ax[1].scatter(
            np.log(self.atomic_transitions["equivalent_width"]/self.atomic_transitions["wavelength"]),
            self.atomic_transitions["abundance"], facecolor="k")

        ax[2].scatter(self.atomic_transitions["wavelength"],
            self.atomic_transitions["profile_stddev"], facecolor="k")

        state, info = self.equalibrium_state(full_output=True)
        coefficients = info["coefficients"]
        ok = np.isfinite(self.atomic_transitions["abundance"])
        x = self.atomic_transitions["excitation_potential"][ok]
        ax[0].plot(x, np.polyval(coefficients[0], x), c='b')
        x = np.log(self.atomic_transitions["equivalent_width"] \
            / self.atomic_transitions["wavelength"])[ok]
        ax[1].plot(x, np.polyval(coefficients[1], x), c='b')
        raise a


    def equalibrium_state(self, data=None, sigma_clip=2, full_output=False):
        """
        Calculate the current equalibrium state.
        """

        if data is None:
            data = self.atomic_transitions

        # Use the abundances in self.atomic_transitions
        finite = np.isfinite(data["abundance"])
        if finite.sum() == 0:
            raise ValueError("no abundances calculated")

        outliers_removed = -1
        neutral = ((data["species"] % 1) == 0)
        ionised = ((data["species"] % 1) > 0)
        metallicity = data.meta["abundances_given_stellar_parameters"][2]

        while True:
            neutral *= finite
            ionised *= finite

            # Calculate the slope with excitation potential and abundance.
            excitation_regression = stats.linregress(
                x=data["excitation_potential"][neutral + ionised],
                y=data["abundance"][neutral + ionised])

            # Calculate the ionisation state.
            ionisation_state = data["abundance"][neutral].mean() \
                - data["abundance"][ionised].mean()

            # Calculate the abundance state.
            abundance_state = np.mean(data["abundance"][neutral + ionised] \
                - (atmospheres.solar_abundance(metallicity \
                    + data["species"][neutral + ionised])))

            # Slope with reduced equivalent width and line abundance.
            rew = np.log(data["equivalent_width"] / data["wavelength"])
            rew_regression = stats.linregress(
                x=rew[neutral + ionised],
                y=data["abundance"][neutral + ionised])

            # Calculate the state
            state = np.array([
                excitation_regression[0],
                ionisation_state,
                abundance_state,
                rew_regression[0]
            ])

            if outliers_removed > -1:
                final_state = (excitation_regression, ionisation_state,
                    abundance_state, rew_regression)
                break

            else:
                # Remove the outliers?
                initial_state = (excitation_regression, ionisation_state,
                    abundance_state, rew_regression)

                if sigma_clip is None or not np.isfinite(sigma_clip) \
                or 0 >= sigma_clip:
                    # Don't remove any outliers
                    outliers_removed, final_state = False, initial_state
                    break

                # We don't want to remove stars that are simply on the edge of
                # distributions because if our initial guess was very wrong, we
                # could just be removing the weakest or strongest lines.

                # Instead we will remove lines that are discrepant from the line
                # fits.

                # Outliers in excitation potential vs abundance:
                x = data["excitation_potential"][neutral + ionised]
                y = data["abundance"][neutral + ionised]
                line = excitation_regression[0] * x + excitation_regression[1]

                differences = np.abs(line - y)
                excitation_sigma = differences/np.std(differences)

                # Outliers in reduced equivalent width vs abundance:
                x = rew[neutral + ionised]
                y = data["abundance"][neutral + ionised]
                line = rew_regression[0] * x + rew_regression[1]

                differences = np.abs(line - y)
                line_strength_sigma = differences/np.std(differences)

                # Update the finite mask to remove outliers
                outliers = (excitation_sigma > sigma_clip) \
                    + (line_strength_sigma > sigma_clip)
                finite[finite] *= ~outliers

                outliers_removed = outliers.sum()
                continue # to re-fit the lines

        # TODO include information about which lines were outliers.
        if full_output:
            info = {
                "initial_state": initial_state,
                "final_state": final_state,
                "coefficients": [
                    [excitation_regression[0], excitation_regression[1]],
                    [rew_regression[0], rew_regression[1]]
                ]
            }
            return (state, info)
        return state


    def estimate_stellar_parameters(self, data, initial_theta=None, **kwargs):
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

        if initial_theta is None:
            initial_theta = self.initial_theta(data)


        # Were profiles measured, and were they measured at the initial theta point?
        profiles_measured = self.atomic_transitions.meta.get(
            "profiles_given_stellar_parameters", False)
        sp_initial_theta = [initial_theta[p] for p in self._stellar_parameters]
        if not profiles_measured or sp_initial_theta != profiles_measured:
            # Measure them again.
            fitted_profiles = self.fit_atomic_transitions(data, **initial_theta)

        # Create a copy of the transitions for our little 'objective adventure'.
        transitions = self.atomic_transitions.copy()
        measured_atomic_lines = np.isfinite(transitions["equivalent_width"]) \
            * (transitions["equivalent_width"] > 0)

        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(2)
        ax[0].scatter(transitions["excitation_potential"], transitions["abundance"],
            facecolor="k")
        ax[1].scatter(np.log(transitions["equivalent_width"]/transitions["wavelength"]),
            transitions["abundance"], facecolor="k")

        def objective_function(theta):

            logger.debug("Attempting theta {}".format(theta))
            stellar_parameters, microturbulence = theta[:3], theta[3]
            
            # Interpolate a photosphere and calculate atomic abundances.
            try:
                photosphere = self._interpolator(*stellar_parameters)
                abundances = synthesis.moog.atomic_abundances(
                    transitions[measured_atomic_lines], photosphere,
                    microturbulence=microturbulence)

            except:
                logger.exception("Exception while calculating abundances at "\
                    "{0}".format(theta))
                # TODO what about when no microturbulence?


                fig, ax = plt.subplots(2)
                ax[0].scatter(transitions["excitation_potential"], transitions["abundance"],
                    facecolor="k")
                ax[1].scatter(np.log(transitions["equivalent_width"]/transitions["wavelength"]),
                    transitions["abundance"], facecolor="k")

                raise WTFError()
                return [np.nan] * 4

            else:
                transitions["abundance"][measured_atomic_lines] = abundances

            state = self.equalibrium_state(data=transitions)
            
            state[2] *= 0.1
            print("theta", theta, state, (state**2).sum())

            # Re-arrange the state.
            return state


        # Optimisation
        xtol = kwargs.pop("xtol", 1e-10)
        maxfev = kwargs.pop("maxfev", 100)

        result = op.fsolve(objective_function, sp_initial_theta,
            fprime=stellar_parameter_jacobian_approximation, col_deriv=1,
            epsfcn=0, xtol=xtol, maxfev=maxfev, full_output=True)

        fig, ax = plt.subplots(2)
        ax[0].scatter(transitions["excitation_potential"], transitions["abundance"],
            facecolor="k")
        ax[1].scatter(np.log(transitions["equivalent_width"]/transitions["wavelength"]),
            transitions["abundance"], facecolor="k")

        raise CaffieneRequiredToContinueError



def wavelengths_in_data(wavelengths, data):
    orders = -np.ones(len(wavelengths), dtype=int)
    for i, spectrum in enumerate(data):
        orders[(spectrum.disp[-1] >= wavelengths) \
            * (wavelengths >= spectrum.disp[0])] = i
    return orders

