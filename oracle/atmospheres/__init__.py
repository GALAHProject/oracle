#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" Model atmospheres """

__author__ = "Andy Casey <arc@ast.cam.ac.uk>"

from .abundances import asplund_2009 as solar_abundance
from .castelli_kurucz import Interpolator as ck_interp
from .marcs import Interpolator as marcs_interp
#from .stagger import Interpolator as stagger_interp
from . import utils


def interpolator(kind=None, **kwargs):

    if kind is None or kind.lower() == "castelli/kurucz":
        return ck_interp(**kwargs)
    elif kind.lower() == "marcs":
        return marcs_interp(**kwargs)