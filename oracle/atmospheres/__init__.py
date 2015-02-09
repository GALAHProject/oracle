#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" Model atmospheres """

__author__ = "Andy Casey <arc@ast.cam.ac.uk>"

from .abundances import asplund_2009 as solar_abundance
# TODO
from .marcs import Interpolator
from . import utils
