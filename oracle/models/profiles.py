# coding: utf-8

""" Absorption profiles """

__author__ = "Andy Casey <arc@ast.cam.ac.uk>"

import numpy as np
from scipy.special import wofz


def gaussian(mu, sigma, x):
    """
    Evaluates a Gaussian profile with ``mu`` and ``sigma`` at all values of ``x``.

    :param mu:
        The profile mean.

    :type mu:
        float

    :param sigma:
        The standard deviation of the profile.

    :type sigma:
        float

    :param x:
        The values to calculate the Gaussian profile at.

    :type x:
        :class:`numpy.array`

    :returns:
        An array with values for the calculated absorption profile.

    :rtype:
        :class:`numpy.array`
    """

    return np.exp(-(x - mu)**2 / (2*sigma**2))


def lorentzian(mu, scale, x):
    """
    Evaluates a Lorentzian absorption profile at all `x` values.

    :param mu:
        The centroid of the line.

    :type mu:
        float

    :param scale:
        The scale parameter (also known as gamma).

    :type scale:
        float

    :param x:
        An array of x-points to calculate the profile at.

    :type x:
        :class:`numpy.array`

    :returns:
        The calculated profile points at ``x``.

    :rtype:
        :class:`numpy.array`
    """

    assert scale >= 0

    # 1./np.pi = 0.3183098861837907
    y = scale/((x - mu)**2 + scale**2)
    return (y - np.min(y))/np.ptp(y)
    

def voigt(mu, fwhm, shape, x):
    """
    Evaluates a Voigt absorption profile across the `x`
    values with a given local continuum.

    V(x,sig,gam) = Re(w(z))/(sig*sqrt(2*pi))
    z = (x+i*gam)/(sig*sqrt(2))

    :param mu:
        The centroid of the line.

    :param fwhm:
        The full-width half-maximum of the Gaussian profile.

    :param shape:
        The shape parameter.

    :param x:
        An array of x-points.
    """
    n = len(x) if not isinstance(x, float) else 1

    profile = 1. / wofz(np.zeros((n)) + 1j * np.sqrt(np.log(2.0)) * shape).real
    return profile * wofz(2*np.sqrt(np.log(2.0)) * (x - mu)/fwhm \
        + 1j * np.sqrt(np.log(2.0))*shape).real