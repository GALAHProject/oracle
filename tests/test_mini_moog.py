#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Make sure the atomic line abundances calculated with compiled MOOG (July 2014) 
version are the same as the hacked Python version.
"""

from __future__ import division, print_function

__author__ = "Andy Casey <arc@ast.cam.ac.uk>"

import numpy as np
import oracle

line_list_filename = "18sco-line-abundances.txt"


def test_18sco(start=0, N=None, debug=True):
    """
    Make sure we get the same abundances as the compiled MOOG version for 18Sco.
    """

    line_list = np.core.records.fromarrays(np.loadtxt(line_list_filename,
        usecols=(0, 1, 2, 3, 4)).T, names=("wavelength", "species", 
        "excitation_potential", "loggf", "equivalent_width"))

    N = len(line_list) if N is None else N
    line_list = line_list[start:start + N]

    # Sort it by wavelength
    indices = np.argsort(line_list["wavelength"])
    line_list = line_list[indices]

    compiled_moog_abundances = np.loadtxt(
        line_list_filename, usecols=(6, ))[start:start + N][indices]

    mini_moog_abundances = oracle.synthesis.moog.atomic_abundances(
        line_list, [5810, 4.44, 0.03], microturbulence=1.07,
        photosphere_kwargs={"kind": "MARCS"}, debug=debug)

    differences = mini_moog_abundances - compiled_moog_abundances
    print("Summary of differences (mean/median/std. dev./|max|): "
        "{0:.2e} / {1:.2e} / {2:.2e} / {3:.2e}".format(
            np.mean(differences), np.median(differences),
            np.std(differences), np.abs(differences).max()))
    assert np.all(np.abs(differences) < 0.005)


def test_repeated_moog_calls(debug=False, rtol=1e-05, atol=1e-08):
    """
    Make the same MOOG call 100 times and make sure that we get the exact
    same abundances each time.
    """

    line_list = np.core.records.fromarrays(np.loadtxt(line_list_filename,
        usecols=(0, 1, 2, 3, 4)).T, names=("wavelength", "species", 
        "excitation_potential", "loggf", "equivalent_width"))

    first_abundances = oracle.synthesis.moog.atomic_abundances(
        line_list, [5810, 4.44, 0.03], microturbulence=1.07,
        photosphere_kwargs={"kind": "MARCS"}, debug=debug)

    for i in range(100):
        nth_abundances = oracle.synthesis.moog.atomic_abundances(
            line_list, [5810, 4.44, 0.03], microturbulence=1.07,
            photosphere_kwargs={"kind": "MARCS"}, debug=debug)

        assert np.allclose(first_abundances, nth_abundances, rtol, atol)
