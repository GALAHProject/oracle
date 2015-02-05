#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" General utilities for parsing Castelli-Kurucz model atmospheres. """

from __future__ import division, absolute_import, print_function

__author__ = "Andy Casey <arc@ast.cam.ac.uk>"

import gzip
import numpy as np


def parse_filename(filename, full_output=False):
    """
    Return the basic stellar parameters from the filename.
    """

    basename = filename.split("/")[-1]
    teff = float(basename.split("t")[1].split("g")[0])
    logg = float(basename.split("g")[1].split("k")[0])/10.
    feh = float(basename[1:4].replace("p", "").replace("m", "-"))/10.
    alpha = [0, 0.4][basename[4] == "a"]
    parameters = [teff, logg, feh, alpha]

    if full_output:
        names = ("effective_temperature", "surface_gravity", "metallicity",
            "alpha_enhancement")
        return (parameters, names)
    return parameters


def parse_photospheric_structure(filename, ndepth=None, line=23,
    full_output=False):
    """
    Parse the photospheric structure (optical depths, temperatures, electron
    and gas pressures) from the filename provided.
    """

    opener = gzip.open if filename[-3:].lower() == ".gz" else open
    with opener(filename, "r") as fp:
        contents = fp.readlines()
    if ndepth is None: # Read it from file if it was not given
        ndepth = int(contents[22].split()[2])

    data = np.array(map(float, 
        "".join(contents[line:line+ndepth]).split())).reshape(ndepth, -1)
    if full_output:
        names = ("RHOX", "T", "P", "XNE", "ABROSS", "ACCRAD", "VTURB", "FLXCNV",
            "VCONV", "VELSND")
        return (data, names)
    return data
