# coding: utf-8

""" Approximate the Jacobian for the equalibria model """

from __future__ import absolute_import, print_function

__all__ = ["approximate", "estimate"]
__author__ = "Andy Casey <arc@ast.cam.ac.uk>"

import logging
import numpy as np

logger = logging.getLogger("oracle")


def approximate(stellar_parameters):
    """
    Calculate the approximate Jacobian matrix, given some stellar parameters.

    :param stellar_parameters:
        The stellar parameters to calculate the Jacobian at. These are expected
        to be the effective temperature, the surface gravity, the overall
        metallicity, and the microturbulence.

    :type stellar_parameters:
        :class:`numpy.array`

    :returns:
        The Jacobian approximation.
    """

    logger.debug("Updating Jacobian at Teff = {0:.0f} K, logg = {1:.2f}, [M/H] "
        "= {2:.2f}, vt = {3:.2f}".format(*stellar_parameters))

    # (effective_temperature, surface_gravity, metallicity, microturbulence)
    m = np.array([
        [ 5.4393e-08,  1.6258e-01,  1.0897e-02, -7.2560e-02],
        [ 4.2613e-08, -5.7948e-01, -1.1533e-01, -4.3985e-01],
        [-3.2710e-08, -1.2006e-01, -2.8592e-05,  3.8185e-03],
        [-1.7822e-08, -1.2114e-01, -1.8847e-02,  3.5564e-02]
    ])
    b = np.array([
        [-4.8623e-04, -8.2654e-01, -2.3837e-02, +1.2853e-01],
        [-4.2039e-04, -1.2402e-00, -9.2341e-02, +8.0592e-02],
        [+2.8178e-04, -3.5816e-02, +1.4257e-03, -1.6601e-02],
        [+1.8250e-04, +4.1779e-01, -1.0949e-01, -1.1024e-01]
    ])

    n = stellar_parameters.size
    x = np.repeat(stellar_parameters, n).reshape(n, n).T
    return (m * x + b).T


def estimate(stellar_parameters, measured_atomic_transitions, spacing=None, 
    atmosphere_kwargs=None, full_output=True):
    """
    Estimate the Jacobian from measured transitions.
    """

    # Sample the full grid at  each spacing.

    # Calculate the state at each point.

    # Approximate the required relations with polynomials of a given order

    # Return an array to dot multiply with the stellar parameters

    raise NotImplementedError