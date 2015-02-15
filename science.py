#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" Analyse the Gaia benchmarks. """

from __future__ import absolute_import, print_function

__author__ = "Andy Casey <arc@ast.cam.ac.uk>"

import os
import logging
import sys
import tarfile
from glob import glob
from time import time
from urllib import urlretrieve

import astropy.table
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import numpy as np

import oracle

logger = logging.getLogger("oracle")

if __name__ != "__main__":
    print("This is a script.")
    sys.exit(0)

# Download the benchmark data and unpack it.
data_url = "https://zenodo.org/record/15103/files/benchmarks.tar.gz"
record_url = "https://zenodo.org/record/15103/"
if not os.path.exists("DATA/benchmarks/benchmarks.csv"):
    logger.info("Downloading {0}".format(data_url))
    try:
        urlretrieve(data_url, "benchmarks.tar.gz")
    except IOError:
        logger.exception("Error downloading benchmark data from {0}".format(data_url))
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

    logger.info("Solving for {0:} ({1:.0f}, {2:.2f}"
            ", {3:.2f})".format(star, benchmark["effective_temperature"],
            benchmark["surface_gravity"], benchmark["metallicity"]))

    # Load the model and the data.
    data = map(oracle.specutils.Spectrum1D.load, filenames)
    model = oracle.models.EqualibriaModel("galah.yaml")

    t_init = time()
    initial_theta = model.initial_theta(data)

    #stellar_parameters = model.estimate_stellar_parameters(data,
    #    initial_theta=initial_theta)

    time_taken = time() - t_init

    results.append([
        star,
        benchmark["effective_temperature"],
        benchmark["surface_gravity"],
        benchmark["metallicity"],
        initial_theta["effective_temperature"],
        initial_theta["surface_gravity"],
        initial_theta["metallicity"],
        time_taken
    ])# + list(stellar_parameters)

    logger.info("Took {0:.1f} seconds to solve for {1}".format(time_taken, star))
    logger.info("Literature for {0:} ({1:.0f}, {2:.2f}"
        ", {3:.2f})".format(star, benchmark["effective_temperature"],
        benchmark["surface_gravity"], benchmark["metallicity"]))
    

# Create a results file in markdown
markdown = \
"""
The [Gaia benchmark spectra]({record_url}) were downloaded and analysed using commit {{commit_sha}}. Literature (leftmost) values are from [Jofre et al. (2014)](http://arxiv.org/pdf/1309.1099v2.pdf).

**Initial Parameters**
Initial stellar parameters, velocities, and continuum coefficients were first estimated by cross-correlation against a grid of models:

Star | Teff | logg | [Fe/H] | Teff [ccf] | logg [ccf] | [Fe/H] [ccf] | Time |
:----|:----:|:----:|:------:|:----------:|:----------:|:------------:|:----:|
     |**(K)**|**(cgs)**|    | **(K)**    | **(cgs)**  |              |**(sec)**|
""".format(record_url=record_url)

for row in results:
    markdown += "{0} | {1:.0f} | {2:.3f} | {3:+.3f} | {4:.0f} | {5:.2f} | {6:+.2f}"\
        " | {7:.1f}\n".format(*row)

# Create a results table for easier plottingg 
results_table = astropy.table.Table(rows=results,
    names=["Star", "Teff_lit", "logg_lit", "[Fe/H]_lit", "Teff_ccf", "logg_ccf",
    "[Fe/H]_ccf", "Time"])

# Make a difference plot
fig, ax = plt.subplots(3)
ax[0].scatter(results_table["Teff_lit"], results_table["Teff_ccf"]-results_table["Teff_lit"], facecolor="k")
#ax[0].scatter(results_table["Teff_lit"], results_table["Teff_eq"]-results_table["Teff_lit"], facecolor="k")
ax[0].axhline(0, ls=":", c="#666666")
ax[0].set_xlabel("$T_{\\rm eff}$ (K)")
ax[0].set_ylabel("$\Delta{}T_{\\rm eff}$ (K)")
_ = np.max(np.abs(ax[0].get_ylim()))
ax[0].set_ylim(-_, +_)
ax[0].yaxis.set_major_locator(MaxNLocator(5))

ax[1].scatter(results_table["logg_lit"], results_table["logg_ccf"]-results_table["logg_lit"], facecolor="k")
#ax[1].scatter(results_table["logg_lit"], results_table["logg_eq"]-results_table["logg_lit"], facecolor="k")
ax[1].axhline(0, ls=":", c="#666666")
ax[1].set_xlabel("$\log{g}$")
ax[1].set_ylabel("$\Delta{}\log{g}$ (dex)")
_ = np.max(np.abs(ax[1].get_ylim()))
ax[1].set_ylim(-_, +_)
ax[1].yaxis.set_major_locator(MaxNLocator(5))

ax[2].scatter(results_table["[Fe/H]_lit"], results_table["[Fe/H]_ccf"]-results_table["[Fe/H]_lit"], facecolor="k")
#ax[2].scatter(results_table["[Fe/H]_lit"], results_table["[Fe/H]_eq"]-results_table["[Fe/H]_lit"], facecolor="k")
ax[2].axhline(0, ls=":", c="#666666")
ax[2].set_xlabel("[Fe/H]")
ax[2].set_ylabel("$\Delta{}{\\rm [Fe/H]}$ (dex)")
_ = np.max(np.abs(ax[2].get_ylim()))
ax[2].set_ylim(-_, +_)
ax[2].yaxis.set_major_locator(MaxNLocator(5))

fig.tight_layout()
fig.savefig("benchmarks.png")

# Try to upload the figure to Imgur
try:
    import pyimgur
    imgur = pyimgur.Imgur(os.environ.get("IMGUR_CLIENT_ID", None))
    uploaded_image = imgur.upload_image("benchmarks.png")

except BaseException as e:
    logger.exception("Could not upload benchmarks image to Imgur")
    markdown += "\n\nError uploading figures to Imgur ({0}: {1})".format(
        e.errno, e.strerror)

else:
    markdown += "\n\n![Benchmark results]({})".format(uploaded_image.link)


# Save the markdown to file.
with open("results.md", "w") as fp:
    fp.write(markdown)

