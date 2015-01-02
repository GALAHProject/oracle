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
from scipy import ndimage, stats, optimize as op

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
        surrounding=1.0, outliers=True, full_output=False):
        """
        Fit an absorption profile to the provided spectrum.
        """

        # allow wavelength to move by some amount?

        # initial guess of fwhm?
        # initial guess of line depth?
        # initial guess of continuum if not given a blending spectrum?

        scalar_continuum = True if continuum is None \
            or isinstance(continuum, (int, float)) else False

        if continuum is None:
            data_range = (wavelength - surrounding, wavelength + surrounding)
            data = spectrum.slice(data_range)
            continuum = 1

        else:
            data_range = (continuum.disp[0], continuum.disp[-1])
            data = spectrum.slice(data_range)
            #continuum = np.interp(data.disp, continuum.disp, continuum.flux)

        index = data.disp.searchsorted(wavelength)
        continuum_at_wavelength = continuum if scalar_continuum else continuum.flux[index]

        if initial_depth is None:
            initial_depth = 1. - data.flux[index]/continuum_at_wavelength

        if initial_fwhm is None:
            mid_line_value = continuum_at_wavelength * (1 - 0.5 * initial_depth)

            # Find the wavelengths from the central wavelength where this value
            # exists
            pos_mid_point, neg_mid_point = np.nan, np.nan
            pos_indices = data.flux[index:] > mid_line_value
            if pos_indices.any():
                pos_mid_point = data.disp[index:][pos_indices][0]

            neg_indices = data.flux[:index] > mid_line_value
            if neg_indices.any():
                neg_mid_point = data.disp[:index][neg_indices][-1]

            # If both positive and negative mid points are nans, then we cannot
            # estimate the FWHM with this method
            if not np.isfinite([pos_mid_point, neg_mid_point]).any():
                raise ValueError("cannot estimate initial fwhm")

            # If either mid point is a nan, take the initial fwhm from one side
            # of the profile
            initial_fwhm = pos_mid_point - neg_mid_point

            if not np.isfinite(initial_fwhm):
                initial_fwhm = np.array([
                    2 * (pos_mid_point - wavelength),
                    2 * (wavelength - neg_mid_point)
                ])
                finite = np.isfinite(initial_fwhm)
                initial_fwhm = initial_fwhm[finite][0]

        # Gaussian case
        def absorption_profile(wavelength, sigma, depth):
            #return ndimage.gaussian_filter(continuum, sigma/np.diff(data.disp).mean())\
            #    * (1 - depth * profiles.gaussian(wavelength, sigma, data.disp))
            if scalar_continuum:
                c_ = continuum

            else:
                smoothed_continuum = ndimage.gaussian_filter(
                    continuum.flux, sigma/np.diff(continuum.disp).mean())
                c_ = np.interp(data.disp, continuum.disp, smoothed_continuum)

            return c_ * (1 - depth * profiles.gaussian(wavelength, sigma,
                data.disp))

        def integrate_absorption_profile(wavelength, sigma, depth):
            return sigma * depth

        # Should we model the outliers?
        if outliers:
            def negative_log_likelihood(theta):
                wavelength, sigma, depth, outlier_mu, outlier_sigma, P = theta

                if not (1 > P > 0) or 0 > sigma or 0 > outlier_sigma \
                or 0 > wavelength:
                    return np.nan

                model_line = absorption_profile(wavelength, sigma, depth)
                model_background = outlier_mu

                ivariance_line = data.ivariance
                ivariance_background = data.ivariance + outlier_sigma**2

                likelihood_line = -0.5 * ((data.flux - model_line)**2 \
                    * ivariance_line - np.log(ivariance_line))
                likelihood_background = -0.5 * ((data.flux - model_background)**2 \
                    * ivariance_background - np.log(ivariance_background))

                likelihood = np.sum(np.logaddexp(
                    np.log(1 - P) + likelihood_line,
                    np.log(P) + likelihood_background))

                # Return the negative likelihood, since this function will be
                # minimised
                return -likelihood

            initial_sigma = initial_fwhm/2.355
            finite = np.isfinite(data.flux)
            initial_outlier_mu = np.median(data.flux[finite])
            initial_outlier_sigma = 1e-4

            # approximate width of the line/approximate width covered by the 
            # input data
            initial_P = 1 - (2. * 5 * initial_sigma)/np.ptp(data.disp)

            # Limit the initial fraction to be between [0, 1]
            tolerance = 1e-4
            initial_P = np.clip(initial_P, tolerance, 1 - tolerance)

            labels = ("wavelength", "sigma", "depth", "outlier_mu",
                "outlier_sigma", "outlier_fraction")

            # Start with some initial parameters
            initial_theta = np.array([wavelength, initial_sigma, initial_depth,
                initial_outlier_mu, initial_outlier_sigma, initial_P])
            initial_text = ", ".join(["{0} = {1:.2f}".format(label, value) \
                for label, value in zip(labels, initial_theta)])
            logger.debug("Initial theta for absorption profile at {0:.2f} A is "\
                "{1}".format(wavelength, initial_text))

            # Optimise the negative log likelihood
            optimal_theta, fopt, num_iter, num_funcalls, warnflag = op.fmin(
                negative_log_likelihood, initial_theta, disp=False,
                full_output=True)

            # Report on the optimised parameters
            optimal_text = ", ".join(["{0} = {1:.2f}".format(label, value) \
                for label, value in zip(labels, optimal_theta)])
            logger.debug("Optimal theta for absorption profile at {0:.2f} A is "\
                "{1}".format(wavelength, optimal_text))

        else:
            # The chi-sq function
            def chi_sq(theta):
                wavelength, sigma, depth = theta
                model = absorption_profile(wavelength, sigma, depth)
                difference = (data.flux - model)**2 * data.ivariance
                return difference[np.isfinite(difference)].sum()

            labels = ("wavelength", "sigma", "depth")

            # Prepare the initial theta values
            initial_sigma = initial_fwhm/2.355
            initial_theta = np.array([wavelength, initial_sigma, initial_depth])

            # Optimise the chi-squared value
            optimal_theta, fopt, num_iter, num_funcalls, warnflag = op.fmin(
                chi_sq, initial_theta, disp=False, full_output=True)

        # Integrate the profile to get an equivalent width
        equivalent_width = integrate_absorption_profile(*optimal_theta[:3]) * 1000. # milliAngstroms
        # [TODO] astropy.units

        # Return result if no additional information is requested
        if not full_output:
            return (optimal_theta, equivalent_width)

        # Get our returning information together
        initial = dict(zip(labels, initial_theta))
        optimal = dict(zip(labels, optimal_theta))

        # Create a smoothed continuum spectrum using the optimal theta
        initially_smoothed_continuum = ndimage.gaussian_filter(continuum.flux,
            initial["sigma"]/np.diff(continuum.disp).mean())
        optimally_smoothed_continuum = ndimage.gaussian_filter(continuum.flux,
            optimal["sigma"]/np.diff(continuum.disp).mean())
        initial_continuum = specutils.Spectrum1D(disp=continuum.disp,
            flux=initially_smoothed_continuum,
            variance=np.zeros(continuum.disp.size))
        optimal_continuum = specutils.Spectrum1D(disp=continuum.disp,
            flux=optimally_smoothed_continuum,
            variance=np.zeros(continuum.disp.size))

        minimised_quantity = "negative_log_likelihood" if outliers else "chi_sq"
        initial_profile = specutils.Spectrum1D(disp=data.disp,
            flux=absorption_profile(wavelength, initial_sigma, initial_depth),
            variance=np.zeros(data.disp.size))
        optimal_profile = specutils.Spectrum1D(disp=data.disp,
            flux=absorption_profile(*optimal_theta[:3]),
            variance=np.zeros(data.disp.size))

        info = {
            "initial_theta": initial,
            "optimal_theta": optimal,
            minimised_quantity: fopt,
            "pixels": data.disp.size,
            "num_iter": num_iter,
            "num_funcalls": num_funcalls,
            "warnflag": warnflag,
            "data": data,
            "initial_continuum": initial_continuum,
            "optimal_continuum": optimal_continuum,
            "initial_profile": initial_profile,
            "optimal_profile": optimal_profile
        }

        return (optimal_theta, equivalent_width, info)



    