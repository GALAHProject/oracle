#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" General utilities for parsing MARCS model atmospheres. """

from __future__ import division, absolute_import, print_function

__author__ = "Andy Casey <arc@ast.cam.ac.uk>"

import gzip
import numpy as np

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
