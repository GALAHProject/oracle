#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" Functions for dealing with MARCS model atmospheres. """

from __future__ import division, absolute_import, print_function

__author__ = "Andy Casey <arc@ast.cam.ac.uk>"

# Standard library.
import gzip
import logging

# Third party.
import numpy as np

# Module-specific.
from .interpolator import Interpolator

# Create logger.
logger = logging.getLogger(__name__)

class MARCSInterpolator(Interpolator):

    def __init__(self):
        return super(self.__class__, self).__init__("marcs-2011-standard.pkl")

    def neighbours(self, *point):
        """
        Return the indices of the neighbouring model points.

        This function will switch between spherical and plane-parallel
        photospheres depending on which photospheres are more available, given
        the stellar parameters required.
        """

        raise NotImplementedError


def parse_filename(filename, full_output=False):
    """
    Return the basic stellar parameters from the filename.
    """

    basename = filename.split("/")[-1]
    teff = basename[1:5]
    logg = basename.split("_")[1][1:]
    feh = basename.split("_")[5][1:]
    parameters = map(float, [teff, logg, feh])

    if full_output:
        names = ("effective_temperature", "surface_gravity", "metallicity")
        return (parameters, names)
    return parameters


def parse_photospheric_structure(filename, ndepth=56, line=25,
    full_output=False):
    """
    Parse the photospheric structure (optical depths, temperatures, electron
    and gas pressures) from the filename provided.
    """

    opener = gzip.open if filename[-3:].lower() == ".gz" else open
    with opener(filename, "r") as fp:
        contents = fp.readlines()

    data = np.array(map(float, 
        "".join(contents[line:line+ndepth]).split())).reshape(ndepth, -1)
    if full_output:
        names = ("logTau5000", "Depth", "T", "Pe", "Pg")
        names = ("k", "lgTauR", "lgTau5", "Depth", "T", "Pe", "Pg", "Prad",
            "Pturb")
        return (data, names)
    return data
