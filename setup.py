# coding: utf-8

""" oracle, the suppository of wisdom """ 

from __future__ import division, absolute_import, print_function

__author__ = "Andy Casey <arc@ast.cam.ac.uk>"

# Standard library.
import os
import re
import subprocess
import sys
from glob import glob
from urllib import urlretrieve

# Third-party.
from numpy.distutils.core import Command, Extension, setup

#try:
#    from setuptools import setup

#except ImportError:
#    from distutils.core import setup

major, minor1, minor2, release, serial =  sys.version_info
open_kwargs = {"encoding": "utf-8"} if major >= 3 else {}

def readfile(filename):
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
        'MyAbfind.f', 'Partfn.f', 'Sunder.f', 'Eqlib.f', 'Nearly.f','Discov.f',
        'Invert.f', 'Gammabark.f', 'Damping.f', 'Lineinfo.f', 'Opacit.f',
        'Blankstring.f', 'Prinfo.f', 'Opacmetals.f', 'Synspec.f', 'Cdcalc.f',
        'Linlimit.f', 'Taukap.f', 'Jexpint.f', 'Partnew.f', 'Opacscat.f',
        'OpacHelium.f', 'OpacHydrogen.f', 'Opaccouls.f', 'Rinteg.f', 
        'Trudamp.f', 'Ucalc.f', 'Voigt.f', 'Fakeline.f', 'Curve.f',
        'Lineabund.f', 'Molquery.f', 'Oneline.f', 
        'Inmodel.f', 'Inlines.f', 'Batom.f', 'Bmolec.f', 'MySynth.f']])

# External data.
if "--with-atmospheres" in map(str.lower, sys.argv):
    atmosphere_paths = [
        ("https://zenodo.org/record/14964/files/castelli-kurucz-2004.pkl",
            "oracle/atmospheres/castelli-kurucz-2004.pkl"),
        ("https://zenodo.org/record/14964/files/marcs-2011-standard.pkl",
            "oracle/atmospheres/marcs-2011-standard.pkl")
    ]
    for url, filename in atmosphere_paths:
        print("Downloading {0} to {1}".format(url, filename))
        try:
            urlretrieve(url, filename)
        except IOError:
            raise("Error downloading file {} -- consider trying without the "
                "--with-atmospheres flag".format(url))
    sys.argv.remove("--with-atmospheres")


setup(
    name="oracle",
    version=version,
    author="Andrew R. Casey",
    author_email="arc@ast.cam.ac.uk",
    packages=["oracle", "oracle.atmospheres", "oracle.models",
        "oracle.specutils", "oracle.synthesis"],
    url="http://www.github.com/andycasey/oracle/",
    description="the suppository of all wisdom",
    long_description=readfile(os.path.join(os.path.dirname(__file__), "README.md")),
    install_requires=readfile(
        os.path.join(os.path.dirname(__file__), "requirements.txt")).split("\n"),
    entry_points={
        "console_scripts": ["oracle = oracle.cli:main"]
    },
    ext_modules=[moog],
    #scripts=["oracle/"],
    include_package_data=True,
    package_data={
        "oracle.atmospheres": ["marcs-2011-standard.pkl", "castelli-kurucz-2004.pkl"],
        "oracle.models": ["galah-ambre-grid.pickle"]
    }
)
