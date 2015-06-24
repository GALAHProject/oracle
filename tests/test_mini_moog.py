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


def test_18sco(N=None):

    line_list = np.core.records.fromarrays(np.loadtxt(line_list_filename,
        usecols=(0, 1, 2, 3, 4)).T, names=("wavelength", "species", 
        "excitation_potential", "loggf", "equivalent_width"))[:N]

    compiled_moog_abundances = np.loadtxt(line_list_filename, usecols=(6, ))[:N]

    mini_moog_abundances = oracle.synthesis.moog.atomic_abundances(
        line_list, [5810, 4.44, 0.03], microturbulence=1.07,
        photosphere_kwargs={"kind": "MARCS"}, debug=True)

    differences = mini_moog_abundances - compiled_moog_abundances
    print("Summary of differences (mean/median/std. dev./|max|): "
        "{0:.3f} / {1:.3f} / {2:.3f} / {3:.3f}".format(
            np.mean(differences), np.median(differences),
            np.std(differences), np.abs(differences).max()))
    assert np.all(np.abs(differences) < 0.001)