#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" oracle, the suppository of all wisdom """ 

from __future__ import division, absolute_import, print_function

__author__ = "Andy Casey <arc@ast.cam.ac.uk>"

# Standard library.
import os
import re
import sys
from urllib import urlretrieve

# Third-party.
from numpy.distutils.core import Extension, setup

major, minor1, minor2, release, serial =  sys.version_info

def readfile(filename):
    open_kwargs = {"encoding": "utf-8"} if major >= 3 else {}
    with open(filename, **open_kwargs) as fp:
        contents = fp.read()
    return contents

version_regex = re.compile("__version__ = \"(.*?)\"")
contents = readfile(os.path.join(
    os.path.dirname(os.path.abspath(__file__)),
    "oracle", "__init__.py"))

version = version_regex.findall(contents)[0]

# Extensions.
moog = Extension(name = "oracle.synthesis._mini_moog",
    sources = ["oracle/synthesis/source/moog/{}".format(each) for each in [
        "MyAbfind.f", "Partfn.f", "Sunder.f", "Eqlib.f", "Nearly.f","Discov.f",
        "Invert.f", "Gammabark.f", "Damping.f", "Lineinfo.f", "Opacit.f",
        "Blankstring.f", "Prinfo.f", "Opacmetals.f", "Synspec.f", "Cdcalc.f",
        "Linlimit.f", "Taukap.f", "Jexpint.f", "Partnew.f", "Opacscat.f",
        "OpacHelium.f", "OpacHydrogen.f", "Opaccouls.f", "Rinteg.f", 
        "Trudamp.f", "Ucalc.f", "Voigt.f", "Fakeline.f", "Curve.f",
        "Lineabund.f", "Molquery.f", "Oneline.f", 
        "Inmodel.f", "Inlines.f", "Batom.f", "Bmolec.f", "MySynth.f"]])

# External data.
if "--with-models" in map(str.lower, sys.argv):
    data_paths = [
        # Model photospheres (Castelli & Kurucz, MARCS, Stagger)
        ("https://zenodo.org/record/14964/files/castelli-kurucz-2004.pkl",
            "oracle/atmospheres/castelli-kurucz-2004.pkl"),
        ("https://zenodo.org/record/14964/files/marcs-2011-standard.pkl",
            "oracle/atmospheres/marcs-2011-standard.pkl"),
        """
        (,
            "oracle/atmospheres/stagger-2013-optical.pkl"),
        (,
            "oracle/atmospheres/stagger-2013-mass-density.pkl"),
        (,
            "oracle/atmospheres/stagger-2013-rosseland.pkl"),
        (,
            "oracle/atmospheres/stagger-2013-height.pkl"),
        """
        # Model spectra (AMBRE grid for GALAH)
        ("https://zenodo.org/record/14977/files/galah-ambre-grid.pkl",
            "oracle/models/galah-ambre-grid.pkl")
    ]
    for url, filename in data_paths:
        print("Downloading {0} to {1}".format(url, filename))
        try:
            urlretrieve(url, filename)
        except IOError:
            raise("Error downloading file {} -- consider trying without the "
                "--with-models flag".format(url))
    sys.argv.remove("--with-models")

# Now the magic.
setup(
    name="oracle",
    version=version,
    author="Andrew R. Casey",
    author_email="arc@ast.cam.ac.uk",
    packages=[
        "oracle", "oracle.atmospheres", "oracle.models", "oracle.specutils",
        "oracle.synthesis"
    ],
    url="http://www.github.com/andycasey/oracle/",
    description="the suppository of all wisdom",
    long_description=readfile(
        os.path.join(os.path.dirname(__file__), "README.md")),
    install_requires=readfile(
        os.path.join(os.path.dirname(__file__), "requirements.txt")).split("\n"),
    entry_points={
        "console_scripts": ["oracle = oracle.cli:main"]
    },
    ext_modules=[moog],
    include_package_data=True,
    package_data={
        "oracle.atmospheres": [
            "marcs-2011-standard.pkl",
            "castelli-kurucz-2004.pkl",
            "stagger-2013-optical.pkl",
            "stagger-2013-mass-density.pkl",
            "stagger-2013-rosseland.pkl",
            "stagger-2013-height.pkl"
        ],
        "oracle.models": ["galah-ambre-grid.pkl"]
    }
)
