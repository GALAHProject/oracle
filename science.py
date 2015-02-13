#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" Analyse the Gaia benchmarks. """

from __future__ import absolute_import, print_function

__author__ = "Andy Casey <arc@ast.cam.ac.uk>"

import os
import logging
import tarfile
from glob import glob
from time import time
from urllib import urlretrieve

import astropy.table
import matplotlib.pyplot as plt
import numpy as np

import oracle

logger = logging.getLogger("oracle")


# Download the benchmark data and unpack it.
if not os.path.exists("DATA/benchmarks/benchmarks.csv"):
    data_uri = "https://zenodo.org/record/15103/files/benchmarks.tar.gz"
    logger.info("Downloading {0}".format(data_uri))
    try:
        urlretrieve(data_uri, "benchmarks.tar.gz")
    except IOError:
        logger.exception("Error downloading benchmark data from {0}".format(data_uri))
        raise
    else:
        with tarfile.open("benchmarks.tar.gz") as tar:
            tar.extractall()

# Load the benchmarks.
results = []
benchmarks = astropy.table.Table.read("DATA/benchmarks/benchmarks.csv")
for benchmark in benchmarks:

    star = benchmark["star"]
    filenames = glob("DATA/benchmarks/{}/*.txt".format(star))
    if len(filenames) == 0:
        logger.info("Skipping {0} because no data files found".format(star))
        continue

    logger.info(
        "Solving for {0:} ({1:.0f}, {2:.2f}"
            ", {3:.2f})".format(star, benchmark["effective_temperature"],
            benchmark["surface_gravity"], benchmark["metallicity"]))

    data = map(oracle.specutils.Spectrum1D.load, filenames)

    t_init = time()
    model = oracle.models.EqualibriaModel("tests/benchmarks/galah.yaml")
    stellar_parameters = model.estimate_stellar_parameters(data)

    t_taken = time() - t_init

    results.append([
        star, benchmark["effective_temperature"], benchmark["surface_gravity"],
        benchmark["metallicity"], model._initial_theta["effective_temperature"],
        model._initial_theta["surface_gravity"], model._initial_theta["metallicity"],
    ] + list(stellar_parameters) + [t_taken])

    logger.info("Took {0:.1f} seconds to solve for {1}".format(t_taken, star))

results = astropy.table.Table(rows=results,
    names=["Star", "Teff_lit", "logg_lit", "[Fe/H]_lit", "Teff_ccf", "logg_ccf",
        "[Fe/H]_ccf", "Teff_eq", "logg_eq", "[Fe/H]_eq", "xi_eq", "Time"])

results["Teff_lit"].unit = "K"
results["Teff_ccf"].unit = "K"
results["Teff_eq"].unit = "K"
results["Time"].unit = "seconds"


# Make a difference plot
fig, ax = plt.subplots(3)
ax[0].scatter(results["Teff_lit"], results["Teff_ccf"]-results["Teff_lit"], facecolor="r")
ax[0].scatter(results["Teff_lit"], results["Teff_eq"]-results["Teff_lit"], facecolor="k")
ax[0].axhline(0, ls=":", c="#666666")
ax[0].set_xlabel("$T_{\\rm eff}$ (K)")
ax[0].set_ylabel("$\Delta{}T_{\\rm eff}$ (K)")
_ = np.max(np.abs(ax[0].get_ylim()))
ax[0].set_ylim(-_, +_)

ax[1].scatter(results["logg_lit"], results["logg_ccf"]-results["logg_lit"], facecolor="r")
ax[1].scatter(results["logg_lit"], results["logg_eq"]-results["logg_lit"], facecolor="k")
ax[1].axhline(0, ls=":", c="#666666")
ax[1].set_xlabel("$\log{g}$")
ax[1].set_ylabel("$\Delta{}\log{g}$ (dex)")
_ = np.max(np.abs(ax[1].get_ylim()))
ax[1].set_ylim(-_, +_)

ax[2].scatter(results["[Fe/H]_lit"], results["[Fe/H]_ccf"]-results["[Fe/H]_lit"], facecolor="r")
ax[2].scatter(results["[Fe/H]_lit"], results["[Fe/H]_eq"]-results["[Fe/H]_lit"], facecolor="k")
ax[2].axhline(0, ls=":", c="#666666")
ax[2].set_xlabel("[Fe/H]")
ax[2].set_ylabel("$\Delta{}{\\rm [Fe/H]}$ (dex)")
_ = np.max(np.abs(ax[2].get_ylim()))
ax[2].set_ylim(-_, +_)

fig.tight_layout()

results.write("benchmark-results.fits")
fig.savefig("benchmark-results.pdf")



