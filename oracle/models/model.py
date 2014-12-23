# coding: utf-8

""" An abstract model class for stellar spectra """

from __future__ import absolute_import, print_function

__all__ = ["Model"]
__author__ = "Andy Casey <arc@ast.cam.ac.uk>"

import cPickle as pickle
import json
import logging
import os
import yaml
import numpy as np
import warnings
from hashlib import md5
from pkg_resources import resource_stream
from functools import partial
from scipy import stats, optimize as op

logger = logging.getLogger("oracle")

from oracle import utils, specutils
from oracle.models import profiles, validation

# Silence 'Polyfit may be poorly conditioned' messages
warnings.simplefilter("ignore", np.RankWarning)


class Model(object):

    def __init__(self, configuration, validate=True):
        """
        A general class to probabilistically model stellar spectra.

        :param configuration:
            The path to the configuration file that contains settings for the
            class, a dictionary describing the configuration, or a string 
            containing the configuration in a YAML-dump format. Only YAML-style
            formats are accepted.

        :type configuration:
            str or dict
        """

        if isinstance(configuration, dict):
            self.config = utils.update_recursively(self.config, configuration)

        else:
            if not os.path.exists(configuration):
                # Probably a string configuration.
                try:
                    supplied_configuration = yaml.load(configuration)

                except:
                    raise IOError("configuration file does not exist or the"\
                        " YAML string provided does not describe a valid "\
                        "dictionary")
                else:
                    # We expect a dictionary.
                    if not isinstance(configuration, dict):
                        raise IOError("configuration file does not exist or the"\
                            " YAML string provided does not describe a valid "\
                            "dictionary")

                    self.config = utils.update_recursively(self.config,
                        supplied_configuration)
            else:
                with open(configuration, "r") as fp:
                    supplied_configuration = yaml.load(fp)

                self.config = utils.update_recursively(self.config,
                    supplied_configuration)

        if validate:
            return validation.validate_configuration(self.config)


    def __str__(self):
        return unicode(self).encode("utf-8")


    def __unicode__(self):
        return u"{module}.Model class".format(module=self.__module__)


    def __repr__(self):
        return u"<{module}.Model object with hash {hash} at {location}>".format(
            module=self.__module__, hash=self.hash[:10], location=hex(id(self)))


    def _load_grid(self):
        """
        This is a temporary function until I can generalise the grid.
        """

        with resource_stream(__name__, "galah-ambre-grid.pickle") as fp:
            grid_points, grid_dispersion, grid_fluxes, grid_pixels = pickle.load(fp)
        grid_fluxes = grid_fluxes.reshape(grid_points.size, sum(grid_pixels))

        return (grid_points, grid_dispersion, grid_fluxes)


    def _continuum_order(self, channel_index):
        """
        Parse the configuration and return the continuum order for some channel
        index.
        """
        
        # If not provided, assume no continuum modelling.
        continuum = self.config["model"].get("continuum", -1)
        if continuum == -1:
            return -1

        # A single value means apply to all observed channels
        if isinstance(continuum, (int, float)):
            return continuum

        # Or it's a list
        try:
            return continuum[channel_index]

        except (KeyError, TypeError):
            logger.warn("Could not interpret continuum order for channel {0} "\
                " ({1}), so returning -1".format(channel_index, continuum))
            return -1


    def _continuum(self, dispersion, channel_index, theta):
        """
        Return the continuum at each dispersion point for the given channel.
        """

        order = self._continuum_order(channel_index)
        if 0 > order:
            return np.ones(dispersion.size)

        coefficients = [theta["continuum.{0}.{1}".format(channel_index, i)] \
            for i in range(order + 1)]

        return np.polyval(coefficients[::-1], dispersion)


    @property
    def hash(self):
        """ Return a MD5 hash of the JSON-dumped model configuration. """ 
        return md5(json.dumps(self.config).encode("utf-8")).hexdigest()


    def mask(self, dispersion, z=0, fill_value=False):
        """
        Return an array mask for a given dispersion array and redshift, based on
        the mask information provided in the model configuration file.

        :param dispersion:
            An array of dispersion values to mask.

        :type dispersion:
            :class:`numpy.array`

        :param z: [optional]
            The redshift to apply.

        :type z:
            float

        :param mask_value: [optional]
            The value to use for the masked dispersion points.

        :type mask_value:
            float-like

        :returns:
            An array of the same length of ``dispersion``, filled with ones,
            except for masked values which contain ``mask_value``.

        :rtype:
            :class:`numpy.array`
        """

        regions = self.config.get("mask", None)
        if regions is not None:
            mask = np.ones(len(dispersion), dtype=bool)
            regions = np.array(regions)
            for start, end in regions * (1. + z):
                indices = dispersion.searchsorted([start, end])
                mask[indices[0]:indices[1]] = fill_value
            return mask

        else:
            return np.ones(len(dispersion), dtype=bool)



    def initial_theta(self, data, full_output=False):
        """
        Return an initial guess of the model parameters theta using no prior
        information.

        :param data:
            The observed data.

        :type data:
            list of :class:`oracle.specutils.Spectrum1D` objects
        """

        # Sort the data from blue to red first
        data = sorted(data, key=lambda x: x.disp[0])

        # Need some cacher object to do fast-on-the-fly-comparisons
        parameters = self.parameters(data)

        # Load the pickled cacher object thing, and work out some estimates
        # of the model parameters
        # Cacher must contain:
        # dict with info about the grid, etc., grid points record array, 
        # dispersion map, gigantic flux array
        grid_points, grid_dispersion, grid_fluxes = self._load_grid()

        theta = {}
        num_pixels = 0
        continuum_coefficients = {}
        expected_channel_disp = []
        expected_channel_fluxes = []
        chi_sqs = np.zeros(grid_points.size)

        # Parallelise channels
        logger.warn("andy you should parallelise this part")
        for i, channel in enumerate(data):

            # Splice the wavelengths
            indices = grid_dispersion.searchsorted([channel.disp[0],
                channel.disp[-1]])

            # Temporarily transform the data to the model dispersion points
            rebinned_channel_disp = grid_dispersion[indices[0]:indices[1]]
            rebinned_channel_flux = np.interp(rebinned_channel_disp,
                channel.disp, channel.flux, left=np.nan, right=np.nan)
            rebinned_channel_ivar = np.interp(rebinned_channel_disp,
                channel.disp, channel.ivariance, left=np.nan, right=np.nan)

            # Reference only finite, unmasked pixels
            useful_pixels = np.isfinite(rebinned_channel_flux) \
                * self.mask(rebinned_channel_disp)

            # And interpolate over the non-useful ones
            if np.any(~useful_pixels):
                rebinned_channel_flux[~useful_pixels] = np.interp(
                    rebinned_channel_disp[~useful_pixels],
                    rebinned_channel_disp[useful_pixels],
                    rebinned_channel_flux[useful_pixels])
                rebinned_channel_ivar[~useful_pixels] = 1e-6

            num_pixels += useful_pixels.sum()

            # Get the continuum order
            order = self._continuum_order(i)

            # Cross-correlate the observed data against the grid
            if ("v_rad" in parameters) or ("v_rad.{}".format(i) in parameters):
                if i == 0:
                    v_rads, v_errs, ccf_peaks = \
                        specutils.cross_correlate.cross_correlate_grid(
                            rebinned_channel_disp,
                            grid_fluxes[:, indices[0]:indices[1]],
                            rebinned_channel_flux.copy(), continuum_order=order,
                            threads=self.config["settings"]["threads"])

                    # Identify one with largest CCF max
                    highest_peak = ccf_peaks.argmax()
                    v_rad = v_rads[ccf_peaks.argmax()]
                    logger.debug("Grid point with highest CCF peak in channel "\
                        "{0} is {1} with v_rad = {2:.1f} km/s (+/- {3:.1f} km/"\
                        "s) and R = {4:.3f}".format(i, utils.readable_dict(
                            grid_points.dtype.names, grid_points[highest_peak]),
                        v_rad, v_errs[highest_peak], ccf_peaks[highest_peak]))

                else:
                    v_rad, v_err, ccf_peak = specutils.cross_correlate.cross_correlate_grid(
                        rebinned_channel_disp, np.array([
                            grid_fluxes[chi_sqs.argmin(), indices[0]:indices[1]]]),
                        rebinned_channel_flux.copy(), continuum_order=order)
                    v_rad, v_err, ccf_peak = v_rad[0], v_err[0], ccf_peak[0]

                    logger.debug("CCF peak in channel {0} is {1} with v_rad = "\
                        "{2:.1f} km/s (+/- {3:.1f} km/s)".format(i, ccf_peak,
                            v_rad, v_err))

                if "v_rad" in parameters:
                    # Global, so add it to a list which we will use to take an
                    # average from later
                    if "v_rad" in theta:
                        theta["v_rad"].append(v_rad)
                    else:
                        theta["v_rad"] = [v_rad]
                else:
                    # Measured radial velocity applies to this channel only
                    theta["v_rad.{}".format(i)] = v_rad

            # Calculate continuum coefficients for each model grid point
            if order >= 0:
                # Note this is assuming you are not doing some Big Ass Matrix(tm)
                # operations.
                dispersion_matrix = np.ones((order+1, rebinned_channel_disp.size))
                for j in range(order + 1):
                    dispersion_matrix[j] *= rebinned_channel_disp**j

                A = (rebinned_channel_flux \
                    / grid_fluxes[:, indices[0]:indices[1]]).T
                coefficients = np.linalg.lstsq(dispersion_matrix.T, A)[0]

                # Save the continuum coefficients (we will use them later)
                continuum_coefficients[i] = coefficients
                continuum = np.dot(coefficients.T, dispersion_matrix)
        
                # Calculate the expected fluxes and the chi-sq value at each point
                expected_fluxes = grid_fluxes[:, indices[0]:indices[1]] * continuum

            else:
                # No continuum treatment.
                expected_fluxes = grid_fluxes[:, indices[0]:indices[1]]

            # Keep the expected fluxes if necessary
            if full_output:
                expected_channel_disp.append(rebinned_channel_disp)
                expected_channel_fluxes.append(expected_fluxes)

            # Add to the chi-sq values
            chi_sqs += np.nansum((rebinned_channel_flux - expected_fluxes)**2\
                    * rebinned_channel_ivar, axis=1)

            # Estimate instrumental broadening
            logger.warn("haven't done instrumental broadening")

            # Determine the best match so far
            grid_index = chi_sqs.argmin()
            r_chi_sq = chi_sqs[grid_index] / (num_pixels - len(parameters) - 1)
            logger.debug(u"Grid point with lowest χ² point so far is {0} with "\
                u"reduced χ² = {1:.1f}".format(utils.readable_dict(
                    grid_points.dtype.names, grid_points[grid_index]), r_chi_sq))

        # Do we need to conglomerate radial velocity measurements together?
        if "v_rad" in parameters:
            median_v_rad = np.median(theta["v_rad"])
            logger.debug("Calculating ensemble radial velocity from {0} to be "\
                "{1} km/s".format(theta["v_rad"], median_v_rad))
            theta["v_rad"] = median_v_rad

        # OK, now determine the nearest point as that which has the lowest chi-sq
        grid_index = chi_sqs.argmin()
        closest_grid_point = grid_points[grid_index]
        r_chi_sq = chi_sqs[grid_index] / (num_pixels - len(parameters) - 1)
        logger.info(u"Grid point with lowest χ² point is {0} with reduced χ² = "\
            "{1:.1f}".format(utils.readable_dict(grid_points.dtype.names,
                closest_grid_point), r_chi_sq))

        # Refine the continuum coefficient estimates using the grid point with
        # the lowest chi-sq value


        # Update theta with the nearest grid point
        theta.update(dict(zip(grid_points.dtype.names, closest_grid_point)))

        # And the continuum coefficients
        for i, coefficients in continuum_coefficients.iteritems():
            for j, coefficient in enumerate(coefficients[:, grid_index]):
                theta["continuum.{0}.{1}".format(i, j)] = coefficient

        # Is microturbulence a parameter?
        if "microturbulence" in parameters:
            theta.setdefault("microturbulence", utils.estimate_microturbulence(
                closest_grid_point["effective_temperature"],
                closest_grid_point["surface_gravity"]))

        # See which parameters could not be estimated
        missing_parameters = set(parameters).difference(theta)
        logger.debug("Could not estimate initial parameters for {}".format(
            ", ".join(missing_parameters)))

        if full_output:
            return (theta, r_chi_sq, np.hstack(expected_channel_disp),
                np.hstack([ecf[grid_index] for ecf in expected_channel_fluxes]))
        return theta


    def fit_absorption_profile(self, wavelength, spectrum, continuum=None,
        initial_fwhm=None, initial_depth=None, wavelength_tolerance=0.10,
        surrounding=1.0):
        """
        Fit an absorption profile to the provided spectrum.
        """

        # allow wavelength to move by some amount?

        # initial guess of fwhm?
        # initial guess of line depth?
        # initial guess of continuum if not given a blending spectrum?

        
        if continuum is None:
            data_range = (wavelength - surrounding, wavelength + surrounding)
            data = spectrum.slice(data_range)
            continuum = 1.

        else:
            data_range = (continuum.disp[0], continuum.disp[-1])
            data = spectrum.slice(data_range)
            continuum = np.interp(data.disp, continuum.disp, continuum.flux)


        index = data.disp.searchsorted(wavelength)

        if initial_depth is None:
            if isinstance(continuum, (float, int)):
                initial_depth = 1. - data.flux[index]/continuum
            else:
                initial_depth = 1. - data.flux[index]/continuum[index]


        if initial_fwhm is None:
            if isinstance(continuum, (float, int)):
                mid_line_value = continuum * (1 - 0.5 * initial_depth)

            else:
                mid_line_value = continuum[index] * (1 - 0.5 * initial_depth)

            # Find the wavelengths from the central wavelength where this value
            # exists
            pos_mid_point = data.disp[index:][(data.flux[index:] > mid_line_value)][0]
            neg_mid_point = data.disp[:index][(data.flux[:index] > mid_line_value)][-1]            
            initial_fwhm = pos_mid_point - neg_mid_point


        def absorption_profile(wavelength, fwhm, depth):
            return continuum * (1 - depth * profiles.gaussian(wavelength, fwhm,
                data.disp))

        import matplotlib.pyplot as plt
        fig, ax = plt.subplots()
        ax.plot(data.disp, data.flux, 'k')
        ax.plot(data.disp, continuum, 'b')
        ax.plot(data.disp, absorption_profile(wavelength, initial_fwhm/2.355, initial_depth), 'r')


        def chi_sq(theta):
            wavelength, sigma, depth = theta
            model = absorption_profile(wavelength, sigma, depth)
            difference = (data.flux - model)**2 * data.ivariance
            return difference[np.isfinite(difference)].sum()

        theta_init = np.array([wavelength, initial_fwhm/2.355, initial_depth])
        theta_opt = op.fmin(chi_sq, theta_init)

        ax.plot(data.disp, absorption_profile(*theta_opt), 'g')

        raise a

        return theta_opt



    