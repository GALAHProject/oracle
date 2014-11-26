# coding: utf-8

from __future__ import absolute_import, print_function

""" Equialibria (excitation & ionisation balance) model for stellar spectra """

__author__ = "Andy Casey <arc@ast.cam.ac.uk>"

import logging
import numpy as np

from oracle import specutils, utils
from oracle.models.model import Model

logger = logging.getLogger("oracle")

class EquialibriaModel(Model):

    """
    A class to forward model stellar spectra. This class performs excitation and
    ionization balances in order to determine stellar parameters.
    """

    # Default configuration for a EqualibriaModel class
    config = {
        "model": {
            "redshift": True,
            "instrumental_resolution": True,
            "continuum": False
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

        super(EquialibriaModel, self).__init__(configuration)
        return None


    def parameters(self, num_data_channels):
        """
        Return the model parameters for some data. The model configuration is
        applicable for a large number of observed channels, so the number of 
        observed channels is required to determine the exact number of model
        parameters.

        :param num_data_channels:
            The number of observed data channels.

        :type num_data_channels:
            int

        :returns:
            A list of model parameters.

        :rtype:
            tuple
        """

        try:
            num_data_channels = int(num_data_channels)
        except (ValueError, TypeError):
            raise TypeError("number of data channels must be an integer")

        if 1 > num_data_channels:
            raise ValueError("number of data channels must be a positive integer")

        parameters = ["effective_temperature", "surface_gravity", "[M/H]",
            "microturbulence"]

        # Single radial velocity for all channels
        if self.config["model"]["redshift"] == True:
            parameters.append("v_rad")

        # Different radial velocity for each channel?
        elif isinstance(self.config["model"]["redshift"], (tuple, list, )):
            parameters.extend(["v_rad.{}".format(i) for i, channel_v_rad in \
                zip(range(num_data_channels), self.config["model"]["redshift"])\
                if channel_v_rad])

        # Instrumental broadening
        if self.config["model"]["instrumental_resolution"]:
            parameters.extend(["instrumental_resolution.{}".format(i) \
                for i in range(num_data_channels)])

        # Continuum treatment
        if isinstance(self.config["model"]["continuum"], (tuple, list)):
            # List contains order for each channel
            for i, order in \
                zip(range(num_data_channels), self.config["model"]["continuum"]):
                parameters.extend(["continuum.{0}.{1}".format(i, j) \
                    for j in range(order + 1)])

        # Atomic line abundances
        atomic_lines = self.config["model"].get("atomic_lines", [])
        if 2 > len(atomic_lines):
            raise ValueError("less than two atomic lines provided for equalibria")
        parameters.extend([self._line_parameter_format.format(*atomic_line[:2:-1]) \
            for atomic_line in atomic_lines])

        return tuple(parameters)


    def initial_theta(self, data, full_output=False):
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
        theta, r_chi_sq, model_dispersion, model_fluxes = super(GenerativeModel,
            self).initial_theta(data, True)

        # Synthesise spectra surrounding all the lines that we have to



        # Now use the continuum parameters 

        missing_parameters = set(self.parameters(len(data))).difference(theta)
        if len(missing_parameters) > 0:
            logger.warn("Missing parameters: {}".format(", ".join(
                missing_parameters)))

        if full_output:
            return (theta, r_chi_sq, model_dispersion, model_fluxes)
        return theta


    def fit(self, data, initial_theta=None, grid=None):

        # Use a grid to do interpolation if provided, otherwise we will have to
        # do radiative transfer on the fly


        return None