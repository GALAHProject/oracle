#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" Absorption profile shapes. """

from __future__ import absolute_import, print_function

__author__ = "Andy Casey <arc@ast.cam.ac.uk>"

__all__ = ["gaussian", "voigt"]

import numpy as np

from scipy import integrate
from scipy.special import wofz


class gaussian(object):

    @staticmethod
    def __call__(x, mu, sigma, amplitude, continuum=1.0):
        """Evaluates a Gaussian absorption profile across the `x` values with a 
        given local continuum.

        continuum : `np.array` or call-able function
            The continuum over the given `x` region.
        """

        try:
            # Is the continuum call-able?
            continuum = continuum(x)
        except TypeError:
            None

        g = amplitude * np.exp(-(x - mu)**2 / (2.0 * sigma**2))
        return continuum - g


    @staticmethod
    def integrate(mu, sigma, amplitude, **kwargs):
        return amplitude * sigma * np.sqrt(2*np.pi)


class voigt(object):

    @staticmethod
    def __call__(x, mu, sigma, amplitude, gamma, continuum=1.0):
        """Evaluates a Voigt absorption profile across the `x`
        values with a given local continuum.

        V(x,sig,gam) = Re(w(z))/(sig*sqrt(2*pi))
        z = (x+i*gam)/(sig*sqrt(2))

        continuum : `np.array` or callable function
            The continuum over the given region.
        """

        try:
            # Is the continuum call-able?
            continuum = continuum(x)
        except TypeError:
            None

        z = ((x - mu) + 1j*gamma) / (sigma * np.sqrt(2))
        v = amplitude * wofz(z).real
        return continuum - v


    @staticmethod
    def integrate(mu, sigma, amplitude, gamma, continuum=1.0, threshold=10):

        f = voigt()
        integral, error = integrate.quad(f, mu - threshold * sigma,
            mu + threshold * sigma, args=(mu, sigma, amplitude, gamma))
        return integral
