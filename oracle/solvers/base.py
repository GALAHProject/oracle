#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" Base solvers to fit spectra. """

from __future__ import absolute_import, print_function

__author__ = "Andy Casey <arc@ast.cam.ac.uk>"

import logging
import numpy as np

logger = logging.getLogger("oracle")

class BaseFitter(object):

    _initialised = False

    def __init__(self, *args, **kwargs):
        self._initialised = True


    @classmethod
    def mask_data(cls, x, y, mask_regions, mask_non_finites=True):
        """
        Return a mask array for the data (x, y).

        Like numpy masked arrays, when the mask is True, the element is said to
        be masked (invalid).
        """

        assert x.size == y.size
        if mask_non_finites:
            mask = ~np.isfinite(y)
        else:
            mask = ~np.ones(x.size, dtype=bool)

        for lower, upper in mask_regions:
            upper = +np.inf if upper is None else upper
            lower = -np.inf if lower is None else lower
            assert upper > lower
            mask += (upper > x) * (x > lower)

        return mask
