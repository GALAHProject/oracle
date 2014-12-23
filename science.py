# coding: utf-8

"""
Science benchmark tests for the Oracle code
-------------------------------------------

Travis will run unit tests after every commit pushed to GitHub. If this commit
is associated with a pull request (e.g., develop -> master) then the .travis.yml
after_install script will run this script.

This script will download all of the required data from Zenodo and analyse them.
If the required environment variables exist, then it will upload the resulting
figures to Amazon Web Services. After the science script has been run, another
script will comment on the GitHub pull request and include the result figures.

"""

from __future__ import print_function

# Download the data we require

# Run the analysis on the benchmark stars
print("science.py done")
