# coding: utf-8

""" A Pythonic interface to MOOG """

from __future__ import absolute_import, print_function

__all__ = ["abundances", "synthesise", "_synthesise"]
__author__ = "Andy Casey <arc@ast.cam.ac.uk>"

import numpy as np

from astropy.table import Table

import oracle.atmospheres
from . import _mini_moog as moog


def _format_transitions(transitions):
    """
    Format input transitions ready for Fortran.

    :param transitions:
        The input transitions table for MOOG. This must contain the wavelength,
        species, excitation potential, loggf, equivalent width, and two optional
        damping coefficients.

    :type transitions:
       :class:`astropy.table.Table` or :class:`numpy.core.recordarray`
    """

    data = transitions if hasattr(transitions, "view") else transitions._data

    columns = ("wavelength", "species", "excitation_potential", "loggf",
        "C6???", "C4???", "equivalent_width")

    transitions_arr = np.zeros((len(data), 7))
    for i, column in enumerate(columns):
        if column not in data.dtype.names: continue
        transitions_arr[:, i] = data[column]
    return np.asfortranarray(transitions_arr)


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


def _format_photosphere(photosphere):

    # Translate the human-readable model photosphere into MOOG-speak.
    moog_modtype = {
        "marcs": "WEBMARCS",
        "castelli/kurucz": "KURUCZ"
    }[photosphere.meta["kind"]]

    d = photosphere if hasattr(photosphere, "view") else photosphere._data
    photosphere_as_arr = np.asfortranarray(d.view(float).reshape(d.size, -1))
    return (moog_modtype, photosphere_as_arr)


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



def atomic_abundances(transitions, photosphere_information, microturbulence=None,
    photospheric_abundances=None, atmosphere_kwargs=None, **kwargs):
    """
    Calculate atomic abundances from measured equivalent widths.

    :param transitions:
        A table containing atomic data for all transitions.

    :type transitions:
        :class:`astropy.table.Table`

    :param photosphere_information:
        This can be a model photosphere or a set of stellar parameters. If a set
        of stellar parameters (Teff, logg, [M/H]) is provided, then a model
        photosphere will be created and supplementary atmosphere information can
        be provided with the `atmosphere_kwargs` argument.

    :type photosphere_information:
        :class:`astropy.table.Table` (model photosphere) or list of float

    :microturbulence: [optional, sometimes]
        The microturbulence for the model atmosphere, in km/s. Microturbulence
        is a required parameter for 1D models, but is not required for <3D>
        models.

    :type microturbulence:
        float

    :param photospheric_abundances: [optional]
        Abundances of chemical elements in the photosphere.

    :type photospheric_abundances:
        :class:`np.array` (TODO update to astropy table)

    :param atmosphere_kwargs: [optional]
        Arguments to supply to the :class:`oracle.atmospheres.Interpolator`
        class, if the `photosphere_information` is a 3-length list of stellar
        parameters. This is ignored if the `photosphere_information` is a
        pre-interpolated model photosphere.

    :type atmosphere_kwargs:
        dict
    """

    debug = kwargs.pop("debug", False)

    # photosphere_information can be a photosphere or a set of stellar parameters
    if not isinstance(photosphere_information, (Table, tuple, list, np.ndarray))\
    or (isinstance(photosphere_information, (tuple, list, np.ndarray)) \
        and len(photosphere_information) != 3):
        raise TypeError("photosphere_information must be an interpolated "
            "photosphere in astropy.table.Table format, or a 3-length list "
            "containing the effective temperature, surface gravity, and "
            "metallicity")

        # We need to interpolate a photosphere.
        if atmosphere_kwargs is None:
            atmosphere_kwargs = {}

        # Pro-Tip: You can provide your own atmosphere interpolator with the
        #          _interpolator keyword argument.
        interpolator = kwargs.pop("_interpolator",
            oracle.atmospheres.Interpolator(**atmosphere_kwargs))
        photosphere = interpolator.interpolate(photosphere_information)
        
    else:
        photosphere = photosphere_information

    modtype, photosphere_arr = _format_photosphere(photosphere)

    # <3D> models do not require microturbulence.
    if modtype == "STAGGER":
        if microturbulence is not None:
            logger.debug("Ignoring microturbulence ({0:.2f}) for {1} models"\
                .format(microturbulence, modtype))
            microturbulence = 0.
    elif microturbulence is None:
        raise ValueError("microturbulence is required for 1D models")

    # Prepare the transitions table.
    transitions = _format_transitions(transitions)
    
    # Prepare the abundance information
    photospheric_abundances = _format_abundances(photospheric_abundances)
    metallicity = photosphere.meta["stellar_parameters"]["metallicity"]

    # Calculate abundances.
    code, output = moog.abundances(metallicity, microturbulence, photosphere_arr,
        photospheric_abundances, transitions, modtype_=modtype, debug_=debug)

    return output

    

