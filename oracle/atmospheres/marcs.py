#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" Functions for dealing with MARCS model atmospheres. """

from __future__ import division, absolute_import, print_function

__author__ = "Andy Casey <arc@ast.cam.ac.uk>"

# Standard library.
import gzip
import logging
from collections import Counter

# Third party.
import numpy as np

# Module-specific.
from oracle.atmospheres.interpolator import BaseInterpolator

# Create logger.
logger = logging.getLogger(__name__)


class Interpolator(BaseInterpolator):

    def __init__(self):
        """
        A class to interpolate spherical and plane-parallel MARCS model
        atmospheres.

        We use the standard-composition 1 Solar mass models with microturbulence
        of 1 km/s in plane-parallel models and 2 km/s in spherical models.

        """
        return super(self.__class__, self).__init__("marcs-2011-standard.pkl")

    def neighbours(self, *point):
        """
        Return the indices of the neighbouring model points.

        This function will switch between spherical and plane-parallel
        photospheres depending on which photospheres are more available, given
        the stellar parameters required.
        """

        # Point is actually length 3, but the fourth column contains
        # 'is_spherical', e.g. 1 for spherical, 0 for plane-parallel.

        # Ignore the faux value:
        point = point[:3]
        assert len(point) == 3, "Expected 3-length point for stellar parameters"

        # Check both spherical and plane-parallel?
        indices = super(self.__class__, self).neighbours(*point)
        sph_or_pp = self.stellar_parameters["is_spherical?"][indices]
        if len(np.unique(sph_or_pp)) > 1:
            # Take whichever has more points.
            is_spherical = Counter(sph_or_pp).most_common()[0][0]
            logger.debug("Selecting {0} models".format(
                ["plane-parallel", "spherical"][int(is_spherical)]))
            indices *= (self.stellar_parameters["is_spherical?"] == is_spherical)

        if 2**len(point) > indices.sum():
            raise ValueError("nope")
        return indices


    def interpolate(self, *point):
        """ 
        Return the interpolated photospheric quantities on a common opacity
        scale.
        """

        # Put in a faux point value that we will ignore later.
        point = [] + list(point) + [np.nan]
        return super(self.__class__, self).interpolate(*point)



def parse_filename(filename, full_output=False):
    """
    Return the basic stellar parameters from the filename.
    """

    basename = filename.split("/")[-1]
    teff = basename[1:5]
    logg = basename.split("_")[1][1:]
    feh = basename.split("_")[5][1:]
    parameters = map(float, [teff, logg, feh, int(basename[0].lower() == "s")])

    if full_output:
        names = ("effective_temperature", "surface_gravity", "metallicity",
            "is_spherical?")
        return (parameters, names)
    return parameters


def parse_photospheric_structure(filename, ndepth=56, line=25,
    columns=("lgTau5", "Depth", "T", "Pe", "Pg"), full_output=False):
    """
    Parse the photospheric structure (optical depths, temperatures, electron
    and gas pressures) from the filename provided.
    """

    opener = gzip.open if filename[-3:].lower() == ".gz" else open
    with opener(filename, "r") as fp:
        contents = fp.readlines()

    all_columns = ["k", "lgTauR", "lgTau5", "Depth", "T", "Pe", "Pg", "Prad",
        "Pturb"]
    data = np.array(map(float, 
        "".join(contents[line:line+ndepth]).split())).reshape(ndepth, -1)

    # Splice by the columns we want.
    indices = np.array([all_columns.index(c) for c in columns])
    data = data[:, indices]

    if full_output:
        return (data, columns)
    return data
