# coding: utf-8

""" General utilities for parsing MARCS model atmospheres """

from __future__ import division, absolute_import, print_function

import gzip
import numpy as np

def parse_filename(filename):
    """ Return the basic stellar parameters from the filename. """

    basename = filename.split("/")[-1]
    teff = basename[1:5]
    logg = basename.split("_")[1][1:]
    feh = basename.split("_")[5][1:]

    return map(float, [teff, logg, feh])


def parse_photospheric_structure(filename, ndepth=56, line=25):
    """
    Parse the photospheric structure (optical depths, temperatures, electron
    and gas pressures) from the filename provided.
    """

    opener = gzip.open if filename[-3:].lower() == ".gz" else open
    with opener(filename, "r") as fp:
        contents = fp.readlines()
    return np.array(map(float, 
        "".join(contents[line:line+ndepth]).split())).reshape(ndepth, -1)
