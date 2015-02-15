#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" Test the equalibrium model. """

from __future__ import absolute_import, print_function

__author__ = "Andy Casey <arc@ast.cam.ac.uk>"

from time import time

import oracle


data = [
    oracle.specutils.Spectrum1D.load("/Users/arc/codes/oracle-with-siu/oracle/tests/data/benchmarks/18Sco/18Sco_narval_blue_noresample.txt"),
    oracle.specutils.Spectrum1D.load("/Users/arc/codes/oracle-with-siu/oracle/tests/data/benchmarks/18Sco/18Sco_narval_green_noresample.txt"),
    oracle.specutils.Spectrum1D.load("/Users/arc/codes/oracle-with-siu/oracle/tests/data/benchmarks/18Sco/18Sco_narval_red_noresample.txt"),
    oracle.specutils.Spectrum1D.load("/Users/arc/codes/oracle-with-siu/oracle/tests/data/benchmarks/18Sco/18Sco_narval_ir_noresample.txt")
]
model = oracle.models.EqualibriaModel("hermes_classical.yaml")

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
stellar_parameters = model.estimate_stellar_parameters(data)
print("Completed in {0:.1f} seconds".format(time() - t_init))


raise a

