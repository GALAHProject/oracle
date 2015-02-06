#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" Convenience script to pickle a set of model atmospheres. """

from __future__ import division, absolute_import, print_function

__author__ = "Andy Casey <arc@ast.cam.ac.uk>"

import cPickle as pickle
import gzip
import os
import sys
from glob import glob

import numpy as np

import marcs
import castelli_kurucz


def pickle_atmospheres(atmosphere_filenames, kind, meta=None):
    """
    Load all model atmospheres, parse the points and photospheric structures.
    """

    if meta is None:
        meta = {
            "kind": kind,
            "source_directory": os.path.dirname(atmosphere_filenames[0])
        }
    elif not isinstance(meta, dict):
        raise TypeError("meta must be a dictionary or None")

    # Get the names from the first filename
    parsers = {
        "marcs": marcs,
        "castelli/kurucz": castelli_kurucz
    }
    try:
        parser = parsers[kind.lower()]
    except KeyError:
        raise ValueError("don't recognise atmosphere kind '{0}'; available kinds"
            " are {1}".format(kind, ", ".join(parsers.keys())))

    _, parameter_names = parser.parse_filename(atmosphere_filenames[0], True)

    # Get the parameters of all the points
    parameters = np.core.records.fromrecords(
        map(parser.parse_filename, atmosphere_filenames), names=parameter_names)

    # Verify there are no duplicates.
    array_view = parameters.view(float).reshape(parameters.size, -1)
    _ = np.ascontiguousarray(array_view).view(np.dtype((np.void,
        array_view.dtype.itemsize * array_view.shape[1])))
    _, idx = np.unique(_, return_index=True)
    
    if idx.size != parameters.size:
        raise ValueError("{} duplicate stellar parameters found".format(
            parameters.size - idx.size))

    # Now sort the array by the left most columns. Keep track of the indices
    # because we will load the photospheres in this order.
    i = np.argsort(parameters, order=parameter_names)
    parameters = parameters[i]

    _, photosphere_columns = parser.parse_photospheric_structure(
        atmosphere_filenames[0], full_output=True)
    d = np.array([parser.parse_photospheric_structure(atmosphere_filenames[_]) \
        for _ in i])

    return (parameters, d, photosphere_columns, meta)



if __name__ == "__main__":

    # Usage: pickler.py <atmosphere_type> <directory> <pickled_filename>

    import argparse

    parser = argparse.ArgumentParser(description="Pickle atmospheres.")
    parser.add_argument("kind", choices=["marcs", "castelli/kurucz"],
        action="store", help="the type of model atmospheres")
    parser.add_argument("directory", action="store", help="directory containing"
        " the atmosphere files")
    parser.add_argument("pickle_filename", action="store",
        help="the filename to save the pickled atmospheres to")

    args = parser.parse_args()

    # Find the files:
    atmosphere_filenames = glob("{}/*".format(args.directory))
    print("Found {0} files in {1}".format(len(atmosphere_filenames),
        args.directory))
    pickled_data = pickle_atmospheres(atmosphere_filenames, args.kind)

    with open(args.pickle_filename, "wb") as fp:
        pickle.dump(pickled_data, fp, -1)
    print("Pickled {0} atmospheres from {1} to {2}".format(args.kind,
        args.directory, args.pickle_filename))

