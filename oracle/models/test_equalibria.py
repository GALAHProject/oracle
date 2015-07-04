#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" Test the equalibrium model. """

from __future__ import absolute_import, print_function

__author__ = "Andy Casey <arc@ast.cam.ac.uk>"

from time import time

import oracle

# 18Sco 5868 89.000 4.230 0.020 -0.860
data = [
    oracle.specutils.Spectrum1D.load("/Users/arc/codes/oracle-with-siu/oracle/tests/data/benchmarks/18Sco/18Sco_narval_blue_noresample.txt"),
    oracle.specutils.Spectrum1D.load("/Users/arc/codes/oracle-with-siu/oracle/tests/data/benchmarks/18Sco/18Sco_narval_green_noresample.txt"),
    oracle.specutils.Spectrum1D.load("/Users/arc/codes/oracle-with-siu/oracle/tests/data/benchmarks/18Sco/18Sco_narval_red_noresample.txt"),
    oracle.specutils.Spectrum1D.load("/Users/arc/codes/oracle-with-siu/oracle/tests/data/benchmarks/18Sco/18Sco_narval_ir_noresample.txt")
]

# Found: 5.95931867e+03   4.96215965e+00  -1.61229813e-01   1.91952005e+00
"""
data = [
    oracle.specutils.Spectrum1D.load("/Users/arc/codes/oracle-with-siu/oracle/tests/data/benchmarks/Arcturus/Arcturus_narval_blue_noresample.txt"),
    oracle.specutils.Spectrum1D.load("/Users/arc/codes/oracle-with-siu/oracle/tests/data/benchmarks/Arcturus/Arcturus_narval_green_noresample.txt"),
    oracle.specutils.Spectrum1D.load("/Users/arc/codes/oracle-with-siu/oracle/tests/data/benchmarks/Arcturus/Arcturus_narval_red_noresample.txt"),
    oracle.specutils.Spectrum1D.load("/Users/arc/codes/oracle-with-siu/oracle/tests/data/benchmarks/Arcturus/Arcturus_narval_ir_noresample.txt")
]
"""

# Found [  4.41904843e+03   2.27098622e+00  -8.92269149e-01   2.09307689e+00]
# 4247,37.000,1.590,0.040,-0.520


import equalibria
model = equalibria.EqualibriaModel("hermes_classical.yaml")

"""
initial_theta, r_chi_sq, expected_dispersion, expected_flux = model.initial_theta(
    data, full_output=True)

fig, axes = plt.subplots(4)
indices = [-1] + list(np.where(np.diff(expected_dispersion) > 10)[0]) + [None]
for i, ax in enumerate(axes):

    disp = expected_dispersion[indices[i]+1:indices[i+1]]
    flux = expected_flux[indices[i]+1:indices[i+1]]

    ax.plot(disp, flux, "bgrr"[i])
    ax.plot(data[i].disp, data[i].flux, 'k')

raise a
"""

t_init = time()
stellar_parameters = model.estimate_stellar_parameters(data, plot=True)
print("Completed in {0:.1f} seconds".format(time() - t_init))


# Upload data to Zenodo

# Script to download data to pwd

# Table of benchmark values

# For each benchmark, check for spectra

# If spectra exists, do the equalibrium

# Put the results into the benchmark table

# Need plotting / visualisation functions or introspection to see what is being fit properly/bad.



raise a

