#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" Model atmospheres """

__author__ = "Andy Casey <arc@ast.cam.ac.uk>"

import logging

from .abundances import asplund_2009 as solar_abundance
from .castelli_kurucz import Interpolator as ck_interp
from .marcs import Interpolator as marcs_interp
#from .stagger import Interpolator as stagger_interp
from . import utils

logger = logging.getLogger("oracle")

def interpolator(kind="castelli/kurucz", **kwargs):

    logger.debug("Initialising {0} atmosphere interpolator".format(kind.upper()))
    kind = kind.lower()
    if kind == "castelli/kurucz":
        return ck_interp(**kwargs)
    elif kind == "marcs":
        return marcs_interp(**kwargs)