# coding: utf-8

from __future__ import absolute_import, print_function

""" Equialibria (excitation & ionisation balance) model for stellar spectra """

__author__ = "Andy Casey <arc@ast.cam.ac.uk>"

import logging
import numpy as np

from oracle import atmospheres, specutils, synthesis, utils
from oracle.models.model import Model

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

    _line_parameter_format = "A({0:.1f}@{1:.2f})"

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
        num_atomic_lines = 0
        wavelength_ranges = [(each.disp[0], each.disp[-1]) for each in data]

        atomic_lines = self.config["model"].get("atomic_lines", None)
        if atomic_lines is None:
            atomic_lines = np.loadtxt(self.config["model"]["atomic_lines_filename"])

        for atomic_line in atomic_lines:
            wavelength, species = atomic_line[:2]

            # Ensure the wavelength is within the observed region
            for wavelength_start, wavelength_end in wavelength_ranges:
                if wavelength_end >= wavelength and wavelength >= wavelength_start:
                    parameters.append(self._line_parameter_format.format(
                        wavelength, species))
                    num_atomic_lines += 1

        if 2 > num_atomic_lines:
            raise ValueError("less than two atomic lines provided for equalibria")

        return tuple(parameters)


    def solve_stellar_parameters(self, measured_atomic_lines, refit_frequency=1,
        jacobian=None, outlier=None):
        """
        Solve for stellar parameters by excitation and ionisation balance.

        """

        if jacobian is None:
            jacobian = jacobian_default

        elif jacobian == False:
            jacobian = None





        return None


    def _unpack_atomic_line(self, transition):
        synthesise_surrounding = 1.0
        opacity_contribution = 1.0
        damp1, damp2 = 0, 0

        wavelength, species, excitation_potential, loggf = transition[:4]
        if len(transition) > 4:
            damp1 = transition[4]

            if len(transition) > 5:
                damp2 = transition[5]

                if len(transition) > 6:
                    synthesise_surrounding = transition[6]

                    if len(transition) > 7:
                        synthesise_surrounding = transition[7]

        return (wavelength, species, excitation_potential, loggf, damp1, damp2,
            synthesise_surrounding, opacity_contribution)


    def _synthesise_blending_spectrum(self, metallicity, microturbulence,
        photospheric_structure, wavelength_start, wavelength_end, wavelength_step,
        opacity_contributes, blending_line_wavelengths=None):

        if blending_line_wavelengths is None:
            blending_line_wavelengths = np.array([l[0] \
                for l in self.config["model"]["blending_lines"]])

        nearby_blending_line = \
            (wavelength_end >= blending_line_wavelengths) * \
            (blending_line_wavelengths >= wavelength_start)

        blending_dispersion = np.arange(wavelength_start, wavelength_end +
            wavelength_step, wavelength_step)

        if blending_line_wavelengths.size == 0 \
        or nearby_blending_line.sum() == 0:
            return (blending_dispersion, np.ones(blending_dispersion.size))

        npoints = blending_dispersion.size
        abundances = np.asfortranarray(np.array([26.0, 7.50]))
        transitions = np.asfortranarray(self.config["model"]["blending_lines"][:, :6])

        syn_limits = np.asfortranarray([wavelength_start, wavelength_end])
        code, blending_dispersion, blending_flux = synthesis.moog._moog.synthesise(
            metallicity, microturbulence, photospheric_structure, abundances,
            transitions, syn_limits, opacity_contributes, npoints_=npoints)

        return (blending_dispersion, blending_flux)


        # Fit an absorption profile in context of the surround
        #    profile_parameters, equivalent_width = self.fit_profile(
        #        wavelength, data[channel_index], 
        #        blending_spectrum=np.vstack([blending_dispersion, blending_flux]).T,
        #        function=self.config["model"]["profile_function"])

    def fit_profile(self, wavelength, data, blending_spectrum=None,
        function="gaussian"):


        raise NotImplementedError


    def initial_theta(self, data, full_output=False, **kwargs):
        """
        Return an initial guess of the model parameters theta using no prior
        information.

        :param data:
            The observed data.

        :type data:
            list of :class:`oracle.specutils.Spectrum1D` objects
        """

        # Get initial theta from the abstract model class, which will estimate
        # radial velocities, continuum coefficients, and stellar parameters
        theta, r_chi_sq, model_dispersion, model_fluxes = super(EqualibriaModel,
            self).initial_theta(data, True)

        # Now we need to estimate abundances of all of the lines:
        # FOR EACH LINE:
        # (1) Take the initial model theta, interpolate a model atmospheres
        # (2) At each line, see if there are any blending lines that we need to
        #     synthesise.
        # (4) The sub-model parameters will be line parameters
        # (5) Fit the sub-model parameters (the profile) and get an equivalent
        #     width
        
        # THEN: Use the equivalent widths to estimate abundances for all of the
        # lines. 

        # OK, now do it:
        # (1) Take the initial model theta, interpolate a model atmospheres
        atmosphere_kwargs = kwargs.pop("atmosphere_kwargs", {})
        interpolator = atmospheres.Interpolator(**atmosphere_kwargs)

        # Get the stellar parameters (microturbulence not necessary for interp)
        parameter_names = ("effective_temperature", "surface_gravity", "[M/H]")
        stellar_parameters = [theta[p] for p in parameter_names]

        # Interpolate the photospheric quantities
        interpolator_kwargs = kwargs.pop("atmosphere_interpolator_kwargs", {})
        photospheric_structure = interpolator(*stellar_parameters, **interpolator_kwargs)


        # Load the wavelengths of the blending lines, since we will access them
        # a lot.
        blending_line_wavelengths = np.array([l[0] \
            for l in self.config["model"].get("blending_lines", [])])

        for i, atomic_line in enumerate(self.config["model"]["atomic_lines"]):

            wavelength, species, excitation_potential, loggf, damp1, damp2, \
                synthesise_surrounding, opacity_contributes = self._unpack_atomic_line(atomic_line)

            # Blending spectrum will need to be synthesised at a sufficiently
            # high sampling rate to match the data. We also need to know *which*
            # channel the transition is in
            pixel_size, channel_indices = minimum_pixel_sampling(data, wavelength)

            if len(channel_indices) > 1:
                logger.warn("Found transition {2} in {1} channels".format(
                    wavelength, len(channel_indices)))
                raise NotImplementedError("haven't decided what to do about this")

            channel_index = channel_indices[0]

            # (3) Synthesise the blending spectrum and apply the continuum to it.
            blending_dispersion, blending_flux = self._synthesise_blending_spectrum(
                photospheric_structure, theta["[M/H]"], theta["microturbulence"], 
                wavelength - synthesise_surrounding,
                wavelength + synthesise_surrounding,
                pixel_size, opacity_contributes,
                blending_line_wavelengths=blending_line_wavelengths)

            # Apply continuum to the blending spectrum
            blending_flux *= self._continuum(blending_dispersion, channel_index,
                theta)

            # Apply radial velocity to the dispersion
            print("vacuum to air pls")
            c, v = 299792.458, theta.get("v_rad.{}".format(channel_index), 0)
            blending_dispersion *= (1. + v/c)

            # Splice a small portion of the spectrum
            continuum = specutils.Spectrum1D(blending_dispersion, blending_flux)

            # Fit an absorption profile in context of the surround
            profile_parameters, equivalent_width = self.fit_absorption_profile(
                wavelength, data[channel_index], continuum=continuum)


            # Multiply the spectrum by the continuum function

            # Apply radial velocity to the model spectrum

            # Update the initial wavelength with the correction for radial
            # velocity

            # Mask the model spectra as appropriate

            # Fit the profile!

        # Create an array of equaivalent widths

        # Calculate abundances from the integrated equalivent width
        #        abundances = synthesis.moog.abundances(stellar_parameters,
        #    microturbulence, measured_atomic_lines)

        # Add them to the measured atomic lines table and turn it into a named
        # array


        missing_parameters = set(self.parameters(data)).difference(theta)
        if len(missing_parameters) > 0:
            logger.warn("Missing parameters: {}".format(", ".join(
                missing_parameters)))

        if full_output:
            return (theta, r_chi_sq, model_dispersion, model_fluxes)
        return theta




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

    return (pixel_size, channel_indices)


def jacobian_default(effective_temperature, surface_gravity, metallicity,
    microturbulence):
    """
    Calculate the approximate Jacobian matrix, given some stellar parameters.
    """

    # Needs teff, vt, logg, feh arrangement
    raise ValueError("update the jacobian matrix!")
    return np.array([
        [ 4.5143e-08*teff - 4.3018e-04, -6.4264e-04*vt + 2.4581e-02, 
            1.7168e-02*logg - 5.3255e-02,  1.1205e-02*feh - 7.3342e-03],
        [-1.0055e-07*teff + 7.5583e-04,  5.0811e-02*vt - 3.1919e-01,
            -6.7963e-02*logg + 7.3189e-02, -4.1335e-02*feh - 6.0225e-02],
        [-1.9097e-08*teff + 1.8040e-04, -3.8736e-03*vt + 7.6987e-03,
            -6.4754e-03*logg - 2.0095e-02, -4.1837e-03*feh - 4.1084e-03],
        [-7.3958e-09*teff + 1.0175e-04,  6.5783e-03*vt - 3.6509e-02,    
            -9.7692e-03*logg + 3.2322e-02, -1.7391e-02*feh - 1.0502e-01]
    ])

