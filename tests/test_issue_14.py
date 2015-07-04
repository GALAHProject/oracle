#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Ensure that Castelli/Kurucz and MARCS model photospheres give similar results for
a single line.
"""

from __future__ import division, print_function

__author__ = "Andy Casey <arc@ast.cam.ac.uk>"

import numpy as np
from oracle import photospheres, synthesis


def test_for_GH_issue_14():

    line = np.core.records.fromarrays(np.array(
        [[4800.648, 26.0, 4.120, -1.028, 63.87]]).T,
        names=("wavelength", "species", "excitation_potential", "loggf",
            "equivalent_width"))

    m = photospheres.interpolator(kind="marcs")
    ck = photospheres.interpolator(kind="castelli/kurucz")

    stellar_parameters, xi = (5777, 4.445, 0), 1.0
    dwarf_ck = ck.interpolate(*stellar_parameters)
    dwarf_marcs = m.interpolate(*stellar_parameters)

    a_ck = synthesis.moog.atomic_abundances(line, dwarf_ck,
        microturbulence=xi)
    a_m = synthesis.moog.atomic_abundances(line, dwarf_marcs,
        microturbulence=xi)

    print(a_ck, a_m, np.mean(np.abs(a_ck - a_m)))
    assert np.all(np.abs(a_ck - a_m) < 0.05)
