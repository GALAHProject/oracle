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
from hashlib import md5
from functools import partial
from scipy import stats

logger = logging.getLogger("oracle")

from oracle import utils, specutils
from oracle.models import validation

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


    @property
    def hash(self):
        """ Return a MD5 hash of the JSON-dumped model configuration. """ 
        return md5(json.dumps(self.config).encode("utf-8")).hexdigest()


    def mask(self, dispersion, z=0, fill_value=np.nan):
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
            mask = np.ones(len(dispersion))
            regions = np.array(regions)
            for start, end in regions * (1. + z):
                indices = dispersion.searchsorted([start, end])
                mask[indices[0]:indices[1]] = fill_value
            return mask

        else:
            return np.ones(len(dispersion), dtype=bool)



    def initial_theta(self, data):
        """
        Return an initial guess of the model parameters theta using no prior
        information.

        :param data:
            The observed data.

        :type data:
            list of :class:`oracle.specutils.Spectrum1D` objects
        """

        # Need some cacher object to do fast-on-the-fly-comparisons
        parameters = self.parameters(len(data))

        # Load the pickled cacher object thing, and work out some estimates
        # of the model parameters
        # Cacher must contain:
        # dict with info about the grid, etc., grid points record array, 
        # dispersion map, gigantic flux array

        with open("galah-points.pickle", "r") as fp:
            grid_points = pickle.load(fp)

        grid_dispersion = np.memmap("galah-blue-dispersion.memmap", mode="r",
            dtype="float32")
        grid_fluxes = np.memmap("galah-blue-fluxes.memmap", mode="r",
            dtype="float32").reshape(grid_points.size, -1)

    
        #grid_description, grid_points, grid_dispersion, grid_fluxes = \
        #    load_grid(grid_filename)

        theta = {}
        num_pixels = 0
        continuum_coefficients = []
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
            finite = np.isfinite(rebinned_channel_flux) \
                * self.mask(rebinned_channel_disp)
            num_pixels += finite.sum()

            # Get the continuum order
            order = self._continuum_order(i)

            # Cross-correlate the observed data against the grid
            if ("v_rad" in parameters) or ("v_rad.{}".format(i) in parameters):
                v_rads, v_errs, ccf_peaks = specutils.cross_correlate.cross_correlate_grid(
                    channel, grid_dispersion[indices[0]:indices[1]],
                    grid_fluxes[:, indices[0]:indices[1]], continuum_order=order,
                    threads=self.config["settings"]["threads"])

                # Identify one with largest CCF max
                highest_peak = ccf_peaks.argmax()
                v_rad = v_rads[ccf_peaks.argmax()]
                logger.info("Grid point with highest CCF peak in channel {0} is"\
                    " {1} with v_rad = {2:.1f} km/s (+/- {3:.1f} km/s) and R = "\
                    "{4:.3f}".format(i, utils.readable_dict(grid_points.dtype.names,
                        grid_points[highest_peak]), v_rad, v_errs[highest_peak],
                    ccf_peaks[highest_peak]))

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
                dispersion_matrix = np.ones((order + 1, finite.size))
                for j in range(order + 1):
                    dispersion_matrix[j] *= rebinned_channel_disp[finite]**j

                A = (rebinned_channel_flux[finite] \
                    / grid_fluxes[:, indices[0]:indices[1]][:, finite]).T
                coefficients = np.linalg.lstsq(dispersion_matrix.T, A)[0]

                # Save the continuum coefficients (we will use them later)
                continuum_coefficients.append(coefficients)
                continuum = np.dot(coefficients.T, dispersion_matrix)
        
                # Calculate the expected fluxes and the chi-sq value at each point
                expected_fluxes = grid_fluxes[:, indices[0]:indices[1]][:, finite]\
                    * continuum

            else:
                # No continuum treatment.
                expected_fluxes = grid_fluxes[:, indices[0]:indices[1]][:, finite]

            # Add to the chi-sq values
            chi_sqs += np.nansum(
                (rebinned_channel_flux[finite] - expected_fluxes)**2\
                    * rebinned_channel_ivar[finite], axis=1)

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

        # Update theta with the nearest grid point
        theta.update(dict(zip(grid_points.dtype.names, closest_grid_point)))

        # Is microturbulence a parameter?
        if "microturbulence" in parameters:
            theta.setdefault("microturbulence", utils.estimate_microturbulence(
                closest_grid_point["effective_temperature"],
                closest_grid_point["surface_gravity"]))

        # See which parameters could not be estimated
        missing_parameters = set(parameters).difference(theta)
        logger.debug("Could not estimate initial parameters for {}".format(
            ", ".join(missing_parameters)))

        #fig, ax = plt.subplots()
        #ax.plot(rebinned_channel_disp, expected_fluxes[grid_index], 'b')
        #ax.plot(channel.disp, channel.flux, 'k')

        return theta
    