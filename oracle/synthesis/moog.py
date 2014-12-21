# coding: utf-8

""" A Pythonic interface to MOOG """

import numpy as np

import oracle.atmospheres
from oracle.synthesis import __moogsilent__ as moogsilent


def synthesise(effective_temperature, surface_gravity, metallicity,
    microturbulence, transitions, wavelength_start, wavelength_end,
    wavelength_step=0.01, opacity_contributes=1.0,
    photospheric_abundances=None, atmosphere_kwargs=None, debug=False):


    # Interpolate the photospheric structure
    if atmosphere_kwargs is None:
        atmosphere_kwargs = {}

    interpolator = oracle.atmospheres.Interpolator(**atmosphere_kwargs)
    photospheric_structure = interpolator.interpolate(effective_temperature,
        surface_gravity, metallicity)

    # Re-arrange the photospheric structure for moogsilent
    photospheric_structure = np.asfortranarray(photospheric_structure
        .view(float).reshape(photospheric_structure.size, -1))

    if photospheric_abundances is None:
        # By default just give the solar Fe ref. to get started.
        abundances = np.array([26.0, 7.50])
    abundances = np.asfortranarray(abundances).reshape(-1, 2)

    syn_limits = np.asfortranarray([wavelength_start, wavelength_end,
        wavelength_step])
    npoints = (wavelength_end - wavelength_start)/wavelength_step + 1

    code, wavelengths, fluxes = moogsilent.synthesise(metallicity, microturbulence,
        photospheric_structure, abundances, transitions, syn_limits,
        opacity_contributes, npoints_=npoints, debug_=debug)

    return (wavelengths, fluxes)


def abundances(effective_temperature, surface_gravity, metallicity,
    microturbulence, transitions, photospheric_abundances=None,
    atmosphere_kwargs=None, debug=False):

    # Transitions should be an array of minimum size, or a record array to
    # indicate what is what

    # Interpolate the photospheric structure
    if atmosphere_kwargs is None:
        atmosphere_kwargs = {}

    interpolator = oracle.atmospheres.Interpolator(**atmosphere_kwargs)
    photospheric_structure = interpolator.interpolate(effective_temperature,
        surface_gravity, metallicity)

    # Re-arrange the photospheric structure for moogsilent
    photospheric_structure = np.asfortranarray(photospheric_structure
        .view(float).reshape(photospheric_structure.size, -1))

    # Prepare the photospheric abundances
    if photospheric_abundances is None:
        photospheric_abundances = np.asfortranarray(np.array([])).reshape(-1, 2)
    else:
        photospheric_abundances = np.asfortranarray(
            np.atleast_2d(photospheric_abundances))

    code, output = moogsilent.abundances(metallicity, microturbulence,
        photospheric_structure, photospheric_abundances, transitions,
        debug_=debug)

    return output

