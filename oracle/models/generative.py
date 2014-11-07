# coding: utf-8

from __future__ import absolute_import, print_function

""" Generative Model for Stellar Spectra """

__author__ = "Andy Casey <arc@ast.cam.ac.uk>"

import cPickle as pickle
import collections
import logging

import numpy as np
from scipy import optimize as op

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

def load_grid(filename):
    with open(filename, "rb") as fp:
        grid_description, grid_points, grid_dispersion, grid_fluxes = \
            pickle.load(fp)
        
    # Reshape the grid fluxes accordingly
    grid_fluxes = grid_fluxes.reshape(-1, grid_dispersion.size)
    return (grid_description, grid_points, grid_dispersion, grid_fluxes)




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

        parameters = ["effective_temperature", "surface_gravity", "metallicity",
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


    def initial_theta(self, data):
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
        theta = super(GenerativeModel, self).initial_guess(data)

        missing_parameters = set(self.parameters(len(data))).difference(theta)
        if len(missing_parameters) > 0:
            logger.warn("Missing parameters: {}".format(", ".join(
                missing_parameters)))

        return theta


    def fit(self, data, initial_theta=None, grid=None):

        # Use a grid to do interpolation if provided, otherwise we will have to
        # do radiative transfer on the fly


        return None


