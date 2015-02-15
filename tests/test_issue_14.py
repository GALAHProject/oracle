
""" Ensure that Castelli/Kurucz and MARCS model atmospheres give similar results
    for a single line. """

import numpy as np
import oracle
from oracle.atmospheres import marcs, castelli_kurucz


def test_for_GH_issue_14():

    line = np.core.records.fromarrays(np.array(
        [[4800.648, 26.0, 4.120, -1.028, 63.87]]).T,
        names=("wavelength", "species", "excitation_potential", "loggf",
            "equivalent_width"))

    m = marcs.Interpolator()
    ck = castelli_kurucz.Interpolator()


    stellar_parameters, xi = (5777, 4.445, 0), 1.0
    dwarf_ck = ck.interpolate(*stellar_parameters)
    dwarf_marcs = m.interpolate(*stellar_parameters)

    a_ck = oracle.synthesis.moog.atomic_abundances(line, dwarf_ck, microturbulence=xi)
    a_m = oracle.synthesis.moog.atomic_abundances(line, dwarf_marcs, microturbulence=xi)

    assert np.all(np.abs(a_ck - a_m) < 0.05)
