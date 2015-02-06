# coding: utf-8

""" A Pythonic interface to MOOG """

from __future__ import absolute_import, print_function

__all__ = ["atomic_abundances", "synthesise"]
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


def _format_photosphere(photosphere_information, photosphere_kwargs,
    interpolator=None):
    """
    Prepare the input photospheric information for MOOG.
    """

    #photosphere_information can be a photosphere or a set of stellar parameters
    if not isinstance(photosphere_information, (Table, tuple, list, np.ndarray))\
    or (isinstance(photosphere_information, (tuple, list, np.ndarray)) \
        and len(photosphere_information) != 3):
        raise TypeError("photosphere_information must be an interpolated "
            "photosphere in astropy.table.Table format, or a 3-length list "
            "containing the effective temperature, surface gravity, and "
            "metallicity")

        # We need to interpolate a photosphere.
        if interpolator is not None:
            if photosphere_kwargs is None:
                photosphere_kwargs = {}
            interpolator = oracle.atmospheres.Interpolator(**photosphere_kwargs)
        photosphere = interpolator.interpolate(photosphere_information)

    else:
        photosphere = photosphere_information

    metallicity = photosphere.meta["stellar_parameters"]["metallicity"]
    modtype = {
        "marcs": "WEBMARCS",
        "castelli/kurucz": "KURUCZ"
    }[photosphere.meta["kind"]]

    d = photosphere if hasattr(photosphere, "view") else photosphere._data
    photosphere_arr = np.asfortranarray(d.view(float).reshape(d.size, -1))

    return (modtype, photosphere_arr, metallicity)


def synthesise(transitions, photosphere_information, wavelength_region,
    wavelength_step=0.01, microturbulence=None, opacity_contribution=1.0,
    photospheric_abundances=None, photosphere_kwargs=None, **kwargs):
    """
    Calculate a synthetic spectrum using the given transitions and photosphere.

    :param transitions:
        A table containing atomic and molecular data for all transitions.

    :type transitions:
        :class:`astropy.table.Table`

    :param photosphere_information:
        This can be a model photosphere or a set of stellar parameters. If a set
        of stellar parameters (Teff, logg, [M/H]) is provided, then a model
        photosphere will be created and supplementary atmosphere information can
        be provided with the `photosphere_kwargs` argument.

    :type photosphere_information:
        :class:`astropy.table.Table` (model photosphere) or list of float

    :param wavelength_region:
        The start and end wavelength to perform the synthesis in. These values
        are expected to be in Angstroms.

    :type wavelength_region:
        2-length tuple of floats

    :microturbulence: [optional, sometimes]
        The microturbulence for the model atmosphere, in km/s. Microturbulence
        is a required parameter for 1D models, but is not required for <3D>
        models.

    :type microturbulence:
        float

    :param wavelength_step: [optional]
        The spacing between synthesis points in Angstroms. Defaults to 0.01 A.

    :type wavelength_step:
        float

    :param opacity_contribution: [optional]
        The maximum distance (in Angstroms) to where each transition contributes
        to the opacity. This defaults to 1 Angstroms.

    :type opacity_contribution:
        float

    :param photospheric_abundances: [optional]
        Abundances of chemical elements in the photosphere.

    :type photospheric_abundances:
        :class:`np.array` (TODO update to astropy table)

    :param photosphere_kwargs: [optional]
        Arguments to supply to the :class:`oracle.atmospheres.Interpolator`
        class, if the `photosphere_information` is a 3-length list of stellar
        parameters. This is ignored if the `photosphere_information` is a
        pre-interpolated model photosphere.

    :type photosphere_kwargs:
        dict
    """

    debug = kwargs.pop("debug", False)
    modtype, photosphere_arr, metallicity = _format_photosphere(
        photosphere_information, photosphere_kwargs,
        interpolator=kwargs.pop("_interpolator", None))

    # <3D> models do not require microturbulence.
    if modtype == "STAGGER":
        if microturbulence is not None:
            logger.debug("Ignoring microturbulence ({0:.2f}) for {1} models"\
                .format(microturbulence, modtype))
            microturbulence = 0.
    elif microturbulence is None:
        raise ValueError("microturbulence is required for 1D models")

    transitions = _format_transitions(transitions)

    # Prepare the abundance information
    photospheric_abundances = _format_abundances(photospheric_abundances)

    if 0 > wavelength_step:
        raise ValueError("wavelength step must be a positive value")

    synthesis_region = np.asfortranarray(
        [] + list(sorted(wavelength_region)) + [wavelength_step])

    pixels = (synthesis_region[1] - synthesis_region[0])/synthesis_region[2] + 1
    code, wavelengths, fluxes = moog.synthesise(metallicity, microturbulence,
        photosphere_arr, photospheric_abundances, transitions, synthesis_region,
        opacity_contribution, npoints_=pixels, modtype_=modtype, debug_=debug)
    return (wavelengths, fluxes)


def atomic_abundances(transitions, photosphere_information, microturbulence=None,
    photospheric_abundances=None, photosphere_kwargs=None, **kwargs):
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
        be provided with the `photosphere_kwargs` argument.

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

    :param photosphere_kwargs: [optional]
        Arguments to supply to the :class:`oracle.atmospheres.Interpolator`
        class, if the `photosphere_information` is a 3-length list of stellar
        parameters. This is ignored if the `photosphere_information` is a
        pre-interpolated model photosphere.

    :type photosphere_kwargs:
        dict
    """

    debug = kwargs.pop("debug", False)
    modtype, photosphere_arr, metallicity = _format_photosphere(
        photosphere_information, photosphere_kwargs,
        interpolator=kwargs.pop("_interpolator", None))

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

    # Calculate abundances.
    code, output = moog.abundances(metallicity, microturbulence, photosphere_arr,
        photospheric_abundances, transitions, modtype_=modtype, debug_=debug)

    return output

    

