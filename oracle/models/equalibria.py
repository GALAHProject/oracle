# coding: utf-8

from __future__ import absolute_import, print_function

""" Equialibria (excitation & ionisation balance) model for stellar spectra """

__author__ = "Andy Casey <arc@ast.cam.ac.uk>"

import logging
import numpy as np
from scipy import stats, optimize as op

from oracle import atmospheres, specutils, synthesis, utils
from oracle.models.model import Model
from oracle.models.jacobian \
    import approximate as stellar_parameter_jacobian_approximation

logger = logging.getLogger("oracle")

class EqualibriaModel(Model):

    """
    A class to forward model stellar spectra. This class performs excitation and
    ionization balances in order to determine stellar parameters.
    """

    # Default configuration for a EqualibriaModel class
    config = {
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
    _line_parameter_format = "A({species:.1f}@{wavelength:.2f})"


    def __init__(self, configuration):
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
        return None


    def parameters(self, data):
        """
        Return the model parameters for some data. The model configuration is
        applicable for a large number of observed channels, so the number of 
        observed channels is required to determine the exact number of model
        parameters.

        :param data:
            The observed data.

        :type num_data_channels:
            list of :class:`oracle.specutils.Spectrum1D` objects

        :returns:
            A list of model parameters.

        :rtype:
            tuple
        """

        parameters = ["effective_temperature", "surface_gravity", "[M/H]",
            "microturbulence"]

        # Single radial velocity for all channels
        if self.config["model"]["redshift"] == True:
            parameters.append("v_rad")

        # Different radial velocity for each channel?
        elif isinstance(self.config["model"]["redshift"], (tuple, list, )):
            parameters.extend(["v_rad.{}".format(i) for i, channel_v_rad in \
                zip(range(len(data)), self.config["model"]["redshift"])\
                if channel_v_rad])

        # Instrumental broadening
        if self.config["model"]["instrumental_resolution"]:
            parameters.extend(["instrumental_resolution.{}".format(i) \
                for i in range(len(data))])

        # Continuum treatment
        if isinstance(self.config["model"]["continuum"], (tuple, list)):
            # List contains order for each channel
            for i, order in \
                zip(range(len(data)), self.config["model"]["continuum"]):
                parameters.extend(["continuum.{0}.{1}".format(i, j) \
                    for j in range(order + 1)])

        # Atomic line abundances
        num_atomic_transitions = 0
        wavelength_ranges = [(each.disp[0], each.disp[-1]) for each in data]

        atomic_transitions = self.config["model"].get("atomic_transitions", None)
        if atomic_transitions is None:
            # If not specified, try and load simplistic information from a file
            atomic_transitions = np.loadtxt(self.config["model"]["atomic_transitions_filename"])

        for atomic_transition in atomic_transitions:

            if isinstance(atomic_transition, dict):
                wavelength = atomic_transition["wavelength"]
                species = atomic_transition["species"]

            else:
                # Assume list
                wavelength, species = atomic_transition[:2]

            # Ensure the wavelength is within the observed region
            for wavelength_start, wavelength_end in wavelength_ranges:
                if wavelength_end >= wavelength and wavelength >= wavelength_start:
                    parameters.append(self._line_parameter_format.format(
                        wavelength=wavelength, species=species))
                    num_atomic_transitions += 1

        if 2 > num_atomic_transitions:
            raise ValueError("less than two atomic lines provided for equalibria")

        return tuple(parameters)


    def initial_theta(self, data, full_output=False, **kwargs):
        """
        Return an initial guess of the model parameters theta using no prior
        information.

        :param data:
            The observed data.

        :type data:
            list of :class:`oracle.specutils.Spectrum1D` objects
        """

        # Get initial theta from the base model class, which will estimate the
        # radial velocities, continuum coefficients, and stellar parameters
        theta, r_chi_sq, model_dispersion, model_fluxes = super(EqualibriaModel,
            self).initial_theta(data, full_output=True, **kwargs)

        if full_output:
            return (theta, r_chi_sq, model_dispersion, model_fluxes)
        return theta


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

        else:
            logger.warn("assert that all parameters are provided by initial theta")

        # Get the initial stellar parameters
        initial_stellar_parameters = [initial_theta[label] for label in \
            ("effective_temperature", "surface_gravity", "[M/H]", 
                "microturbulence")]

        # Optimisation
        state_kwargs = kwargs.copy()
        state_kwargs.update({
            "theta": initial_theta,
            "full_output": False
        })
        state_function = lambda x: self.equalibrium_state(data, *x,
            **state_kwargs)

        # Optional parameters
        xtol = kwargs.pop("xtol", 1e-10)
        maxfev = kwargs.pop("maxfev", 50)

        # Optimise the state function
        result = op.fsolve(state_function, initial_stellar_parameters,
            fprime=stellar_parameter_jacobian_approximation, col_deriv=True,
            epsfcn=0, xtol=xtol, maxfev=maxfev, full_output=True)



        raise NotImplementedError


    def equalibrium_state(self, data, effective_temperature,
        surface_gravity, metallicity, microturbulence, theta=None,
        full_output=False, **kwargs):
        """
        Return the excitation and ionisation state. This information includes:

        # slope of mean line abundance vs excitation excitation_potential
        # mean difference between neutral and ionised lines
        # slope of mean line abundance vs reduced equivalent width
        # mean difference between input metallicty and mean output metallicity

        """

        if theta is None:
            theta = {}

        interpolator = kwargs.pop("interpolator", None)
        if interpolator is None:
            # Interpolate a model atmospheres
            atmosphere_kwargs = kwargs.pop("atmosphere_kwargs", {})
            interpolator = atmospheres.Interpolator(**atmosphere_kwargs)

        # Interpolate the photospheric quantities
        interpolator_kwargs = kwargs.pop("atmosphere_interpolator_kwargs", {})
        photospheric_structure = interpolator(effective_temperature,
            surface_gravity, metallicity, **interpolator_kwargs)        

        # We want to save the profile information and the details of the
        # measured atomic transitions
        profiles = []
        measured_atomic_transitions = []
        for i, atomic_transition in enumerate(self.config["model"]["atomic_transitions"]):

            wavelength, species, excitation_potential, loggf, van_der_waals_broadening, damp2, \
                synthesise_surrounding, opacity_contributes = utils.unpack_atomic_transition(atomic_transition)

            # Blending spectrum will need to be synthesised at a sufficiently
            # high sampling rate to match the data. We also need to know *which*
            # channel the transition is in

            try:
                pixel_size, channel_indices = minimum_pixel_sampling(
                    data, wavelength)

            except ValueError:
                logger.exception("Cannot find wavelength {0:.3f} in any data "
                    "channel. Skipping this transition.".format(wavelength))
                continue
                

            if len(channel_indices) > 1:
                logger.warn("Found transition {2} in {1} channels".format(
                    wavelength, len(channel_indices)))
                raise NotImplementedError("haven't decided what to do about this")

            channel_index = channel_indices[0]

            # Synthesise the blending spectrum
            blending_dispersion, blending_flux = synthesis._synthesise(
                photospheric_structure, metallicity, microturbulence,
                atomic_transition.get("blending_transitions", None),
                wavelength - synthesise_surrounding,
                wavelength + synthesise_surrounding,
                wavelength_step=pixel_size/4., # Oversample the pixels
                opacity_contributes=opacity_contributes)

            # Apply continuum to the blending spectrum
            blending_flux *= self._continuum(blending_dispersion, channel_index,
                theta)

            # Apply radial velocity to the dispersion
            c, v = 299792.458, theta.get("v_rad.{}".format(channel_index), 0)
            blending_dispersion *= (1. + v/c)

            # [TODO] astropy.constants

            # Create a spectrum object of the continuum
            continuum = specutils.Spectrum1D(blending_dispersion, blending_flux)

            # Fit an absorption profile in context of the surrounding region
            try:
                profile_parameters, equivalent_width, profile_info = \
                    self.fit_absorption_profile(wavelength, data[channel_index],
                        continuum=continuum, full_output=True)

            except ValueError:
                logger.exception("Failed to measure atomic transition at "
                    "{0:.3f}:".format(wavelength))
                continue

            else:
                # Save the information
                profiles.append(profile_info)
                measured_atomic_transitions.append([
                    wavelength, species, excitation_potential, loggf,
                    van_der_waals_broadening, damp2, equivalent_width])

        # Calculate abundances from the integrated equalivent widths
        # [TODO] Use the pre-interpolated photospheric structure?
        measured_atomic_transitions = np.array(measured_atomic_transitions)

        # Only use positive lines
        minimum_equivalent_width = np.clip(
            kwargs.pop("minimum_equivalent_width", 5), 0, np.inf)
        filter_lines = (measured_atomic_transitions[:, 6] > minimum_equivalent_width)

        abundances = synthesis.abundances(effective_temperature,
            surface_gravity, metallicity, microturbulence,
            transitions=measured_atomic_transitions[filter_lines])

        # Add them to the measured atomic lines table and turn it into a recarray
        acceptable_atomic_transitions = np.core.records.fromarrays(
            np.vstack([np.array(measured_atomic_transitions[filter_lines]).T,
                abundances]),
            names=("wavelength", "species", "excitation_potential", "loggf",
                "van_der_waals_broadening", "damp2", "equivalent_width", "abundance"))

        # Calculate the excitation and ionisation state
        logger.warn("Assuming only one kind of atomic species")
        neutral_lines = (acceptable_atomic_transitions["species"] % 1) == 0
        ionised_lines = ~neutral_lines

        outliers_removed = False
        while True:

            # Excitation slope 
            excitation_slope = stats.linregress(
                x=acceptable_atomic_transitions["excitation_potential"][neutral_lines],
                y=acceptable_atomic_transitions["abundance"][neutral_lines])[0]

            # Slope with reduced equivalent width and line abundance
            reduced_equivalent_width = np.log(
                acceptable_atomic_transitions["equivalent_width"]\
                    /acceptable_atomic_transitions["wavelength"])
            line_strength_slope = stats.linregress(
                x=reduced_equivalent_width[neutral_lines],
                y=acceptable_atomic_transitions["abundance"][neutral_lines])[0]

            # Calculate the ionisation state
            ionisation_state = \
                acceptable_atomic_transitions["abundance"][neutral_lines].mean() \
                    - acceptable_atomic_transitions["abundance"][ionised_lines].mean()

            # Calculate the abundance state
            abundance_state = (acceptable_atomic_transitions["abundance"] \
                - (atmospheres.solar_abundances(acceptable_atomic_transitions["species"])
                    + metallicity)).mean()

            if outliers_removed:
                # Re-fit the state
                break

            else:
                # Remove the outliers
                outliers_removed = True

        # Collate the state information together (temperature, surface gravity,
        # metallicity, microturbulence)
        state = np.array([
            excitation_slope,
            ionisation_state,
            abundance_state,
            line_strength_slope
        ])

        import matplotlib.pyplot as plt
        for i, (profile, use) in enumerate(zip(profiles, filter_lines)):
            #if not use: continue
            fig, ax = plt.subplots()
            ax.plot(profile["data"].disp, profile["data"].flux, 'k')
            ax.plot(profile["initial_continuum"].disp, profile["initial_continuum"].flux, 'r:')
            ax.plot(profile["initial_profile"].disp, profile["initial_profile"].flux, 'r')

            ax.plot(profile["optimal_continuum"].disp, profile["optimal_continuum"].flux, 'g:')
            ax.plot(profile["optimal_profile"].disp, profile["optimal_profile"].flux, 'g')
            ax.set_title(str(use))
            ax.set_xlim(profile["data"].disp[0], profile["data"].disp[-1])

            filename = "fig-{0:.2f}.png".format(profile["initial_theta"]["wavelength"])
            fig.savefig(filename)
            plt.close("all")
            print("created {0}".format(filename))

        for filtered, profile in zip(filter_lines, profiles):
            profile["filtered"] = filtered

        raise a

        if full_output:
            return (state, acceptable_atomic_transitions, profiles)
        return state




def minimum_pixel_sampling(data, wavelength):
    """
    Find the minimum pixel size in the data at the given wavelength.
    """

    pixel_size = np.nan
    channel_indices = []
    for i, spectrum in enumerate(data):
        if spectrum.disp[-1] >= wavelength >= spectrum.disp[0]:
            channel_indices.append(i)

            index = spectrum.disp.searchsorted(wavelength)

            diff = np.diff(spectrum.disp)
            if pixel_size > diff[index] or not np.isfinite(pixel_size):
                pixel_size = diff[index]

    if len(channel_indices) == 0 and not np.isfinite(pixel_size):
        raise ValueError("cannot find wavelength {0:.2f} in any data channel"\
            .format(wavelength))

    return (pixel_size, channel_indices)

