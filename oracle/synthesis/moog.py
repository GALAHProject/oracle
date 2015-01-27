# coding: utf-8

""" A Pythonic interface to MOOG """

from __future__ import absolute_import, print_function

__all__ = ["abundances", "synthesise", "_synthesise"]
__author__ = "Andy Casey <arc@ast.cam.ac.uk>"

import numpy as np

import oracle.atmospheres
#from . import _moog as moogsilent
from . import _mini_moog as moog
#import oracle.synthesis._moog as moogsilent


def _format_transitions(transitions):
    """
    Format input transitions ready for Fortran.

    :param transitions:
        A list/array of input transitions. This must contain (at least) the
        rest wavelength, species, excitation potential and oscillator strength.
        It can also include two damping coefficients and the equivalent width.

    :type transitions:
        list, :class:`numpy.ndarray`, or :class:`numpy.core.recordarray`
    """

    transitions = np.atleast_2d(transitions)
    n_transitions, n_columns = transitions.shape

    formatted_transitions = np.zeros((n_transitions, 7))
    formatted_transitions[:, :n_columns] = transitions[:]
    return np.asfortranarray(formatted_transitions)


def _format_abundances(abundances=None):
    """
    Format input phospheric abundances ready for Fortran.

    :param abundances: [optional]
        The input photospheric abundances. If no abundances are supplied, then
        the solar composition is assumed.

    :type abundances:
        :class:`numpy.array`
    """

    if abundances is None or len(abundances) == 0:
        formatted_abundances = np.asfortranarray(np.array([])).reshape(-1, 2)
        
    else:
        formatted_abundances = np.asfortranarray(np.atleast_2d(
            photospheric_abundances))

    return formatted_abundances


def _synthesise(photospheric_structure, metallicity, microturbulence,
    transitions, wavelength_start, wavelength_end, wavelength_step=0.01,
    opacity_contributes=1.0, photospheric_abundances=None, oversample=1,
    debug=False):
    """
    Synthesise a spectrum given some photospheric structure, metallicity,
    microturbulence, and array of transitions.
    """

    # If no transitions are provided, just return a normalised continuum
    if transitions is None or len(transitions) == 0:
        wls = np.arange(wavelength_start, wavelength_end + wavelength_step,
            wavelength_step)
        return (wls, np.ones(wls.size))

    oversample = int(oversample)
    if 1 > oversample:
        raise ValueError("oversampling rate must be greater than 1")

    # Format the arrays as necessary
    photospheric_structure = np.asfortranarray(photospheric_structure
        .view(float).reshape(photospheric_structure.size, -1))
    photospheric_abundances = _format_abundances(photospheric_abundances)
    transitions = _format_transitions(transitions)

    # Prepare synthesis limits
    delta = wavelength_step/oversample
    syn_limits = np.asfortranarray([wavelength_start, wavelength_end, delta])
    npoints = (wavelength_end - wavelength_start)/delta + 1

    code, wavelengths, fluxes = moog.synthesise(
        metallicity, microturbulence, photospheric_structure,
        photospheric_abundances, transitions, syn_limits, opacity_contributes,
        npoints_=npoints, debug_=debug)

    return (wavelengths, fluxes)


def synthesise(effective_temperature, surface_gravity, metallicity,
    microturbulence, transitions, wavelength_start, wavelength_end,
    wavelength_step=0.01, opacity_contributes=1.0, photospheric_abundances=None,
    oversample=1, atmosphere_kwargs=None, debug=False):

    if atmosphere_kwargs is None:
        atmosphere_kwargs = {}

    oversample = int(oversample)
    if 1 > oversample:
        raise ValueError("oversampling rate must be greater than 1")

    # Interpolate the photospheric structure
    interpolator = oracle.atmospheres.Interpolator(**atmosphere_kwargs)
    photospheric_structure = interpolator.interpolate(effective_temperature,
        surface_gravity, metallicity)

    # Format the arrays as necessary
    photospheric_structure = np.asfortranarray(photospheric_structure
        .view(float).reshape(photospheric_structure.size, -1))
    photospheric_abundances = _format_abundances(photospheric_abundances)
    transitions = _format_transitions(transitions)

    # Prepare synthesis limits
    delta = wavelength_step/oversample
    syn_limits = np.asfortranarray([wavelength_start, wavelength_end, delta])
    npoints = (wavelength_end - wavelength_start)/delta + 1

    code, wavelengths, fluxes = moog.synthesise(metallicity, microturbulence,
        photospheric_structure, photospheric_abundances, transitions, syn_limits,
        opacity_contributes, npoints_=npoints, debug_=debug)

    return (wavelengths, fluxes)


def abundances(effective_temperature, surface_gravity, metallicity,
    microturbulence, transitions, photospheric_abundances=None,
    atmosphere_kwargs=None, debug=False):

    if 0 >= effective_temperature:
        raise ValueError("effective temperature must be a positive quantity")

    if 0 > microturbulence:
        raise ValueError("microturbulence must be a positive quantity")

    # Transitions should be an array of minimum size, or a record array to
    # indicate what is what

    if atmosphere_kwargs is None:
        atmosphere_kwargs = {}

    # Interpolate the photospheric structure
    interpolator = oracle.atmospheres.Interpolator(**atmosphere_kwargs)
    photospheric_structure = interpolator.interpolate(effective_temperature,
        surface_gravity, metallicity)

    # Format the arrays as necessary
    photospheric_structure = np.asfortranarray(photospheric_structure
        .view(float).reshape(photospheric_structure.size, -1))
    photospheric_abundances = _format_abundances(photospheric_abundances)
    transitions = _format_transitions(transitions)

    code, output = moog.abundances(metallicity, microturbulence,
        photospheric_structure, photospheric_abundances, transitions,
        debug_=debug)

    return output

