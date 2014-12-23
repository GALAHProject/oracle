# coding: utf-8

from __future__ import absolute_import, print_function

""" Generative Model for Stellar Spectra """

__author__ = "Andy Casey <arc@ast.cam.ac.uk>"

import logging
import numpy as np

from oracle import specutils, utils
from oracle.models.model import Model

logger = logging.getLogger("oracle")

# Process is:
#model = models.GenerativeModel("my_filename.yaml")
# model filename contains:
# masks to use
# normalisation to use (in terms of a list of rules per observed channel)
# redshift to apply: single/different redshifts per channel
# instrumental broadening: yes/no since this is assumed to be different per arm
# radiative transfer code to use
# line list location to use
# model atmospheres to use.
# elemental abundances to consider

# model parameters are determined when data is fed to it, since it depends
# on *how* many channels are provided,
#data = []
#initial_theta, initial_r_chi_sq, info = model.initial_theta(data)
#optimised_theta, model.fit(data, initial_theta=What)

class GenerativeModel(Model):

    """
    A class to forward model stellar spectra. This class does multi-dimensional
    interpolation of model spectra and on-the-fly synthesis to generate stellar
    spectra for some given set of model parameters theta.
    """

    # Default configuration for a GenerativeModel class
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

    def __init__(self, configuration):
        """
        Initialise a GenerativeModel class.
        
        :param configuration:
            The path to the configuration file that contains settings for the
            class, a dictionary describing the configuration, or a string 
            containing the configuration in a YAML-dump format. Only YAML-style
            formats are accepted.

        :type configuration:
            str or dict
        """

        super(GenerativeModel, self).__init__(configuration)
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

        # The Generative model does not strictly require the data, only the 
        # number of observed data. But this is to be consistent with other
        # models.
        num_data_channels = len(data)
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



