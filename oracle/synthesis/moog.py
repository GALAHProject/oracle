# coding: utf-8

""" A Pythonic interface to MOOG """

from __future__ import absolute_import, print_function

__author__ = "Andy Casey <arc@ast.cam.ac.uk>"
__all__ = ["atomic_abundances", "synthesise"]

import logging
import multiprocessing
import numpy as np
import os

from astropy.table import Table

import oracle.atmospheres
from . import _mini_moog as moog

logger = logging.getLogger("oracle")


class MOOGException(BaseException):
    def __call__(self, status="MOOG fell over unexpectedly"):
        raise self.__class__(status)


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

    d = transitions if hasattr(transitions, "view") else transitions.as_array()

    columns = (
        "wavelength",
        "species",
        "excitation_potential",
        "loggf",
        #"C1", # RADIATION DAMPING
        "C6", # Van Der Waals DAMPING
        "D0", # Dissociation energy [eV]
        #"C4", # GAMMA, QUADRATIC STARK DAMPING
        "equivalent_width")

    if np.any((d["species"] % 1) > 0.15):
        # [TODO] Maybe we should actually raise an exception here.
        logger.warn("Species with high ionisation states (>1) detected!")

    transitions_arr = np.zeros((len(d), 7), dtype=float)
    for i, column in enumerate(columns):
        if column not in d.dtype.names: continue
        transitions_arr[:, i] = d[column]

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
        formatted_abundances = np.asfortranarray(np.atleast_2d(abundances))

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

    if isinstance(photosphere_information, (tuple, list, np.ndarray)):
        # We need to interpolate a photosphere.
        if interpolator is None:
            if photosphere_kwargs is None:
                photosphere_kwargs = {}
            interpolator = oracle.atmospheres.interpolator(**photosphere_kwargs)
        photosphere = interpolator.interpolate(*photosphere_information)

    else:
        photosphere = photosphere_information

    d = photosphere if hasattr(photosphere, "view") else photosphere.as_array()
    photosphere_arr = np.asfortranarray(d.view(float).reshape(d.size, -1))

    kind = photosphere.meta["kind"].lower()
    if kind == "marcs":

        # Photospheric quantities and units:
        # log(tau) optical depth
        # T        temperature at {tau} (K)
        # Pe       electron pressure at {tau} [dyne/cm^2]
        # Pg       gas pressure at {tau} [dyne/cm^2]

        # Photospheric quantities expected by MOOG:
        # tauref, t, ne, pgas

        modtype = "WEBMARCS"
        indices = np.array([photosphere.dtype.names.index(c) \
            for c in ("lgTau5", "T", "Pe", "Pg")]) # MOOG calls Pe as Ne??
        photosphere_arr = photosphere_arr[:, indices]
        
    elif kind == "castelli/kurucz":
        
        # Photospheric quantities and units:
        # RHOX  density at {tau}
        # T     Temperature at {tau} (K)
        # P     Gas pressure at {tau} [dyne/cm^2]
        # XNE   Number density of electrons at {tau}??? (???)  
        # kappa Rosseland mean opacity on a mass scale (cm^2/gm)

        # Photospheric quantities expected by MOOG:
        # rhox, t, pgas, ne, kaprefmass

        modtype = "KURUCZ"
        indices = np.array([photosphere.dtype.names.index(c) \
            for c in ("RHOX", "T", "P", "XNE", "ABROSS")])
        photosphere_arr = photosphere_arr[:, indices]

    elif kind == "stagger":

        # Photospheric quantities and units:
        # log(tau) optical depth (at 500 nm)
        # rho      density at {tau} (g/cm^3)
        # T        temperature at {tau} (K)
        # Pth      thermodynamic pressure at {tau} [dyne/cm^2]
        # Ptb      turbulent pressure at {tau} [dyne/cm^2]
        # Pe       electron pressure at {tau} [dyne/cm^2]

        # Note that Ptb is not relevant/required for photospheres averaged along
        # the Rosseland opacity.

        # Photospheric quantities passing to MOOG (as 'WEBMARCS' modtype):
        # tauref, t, ne, pgas

        # Note Stagger Pth ~= MARCS Pg

        averaging = photosphere.meta["horizontal_averaging"].lower()
        if averaging[0] == "r":
            raise NotImplementedError("rosseland opacity averages not set up in"
                " MOOG yet")

        else:
            modtype = "WEBMARCS"
            indices = np.array([photosphere.dtype.names.index(c) \
                for c in ("logtau", "T", "Pe", "Pth")])
            photosphere_arr = photosphere_arr[:, indices]

    else:
        raise ValueError("photosphere kind {} not recognised".format(kind))

    metallicity = photosphere.meta["stellar_parameters"]["metallicity"]
    return (modtype, photosphere_arr, metallicity)


def synthesise(transitions, photosphere_information, wavelength_region=None,
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

    :param wavelength_region: [optional]
        The start and end wavelength to perform the synthesis in. These values
        are expected to be in Angstroms. If not specified, then the region that
        is +/- `opacity_contribution` around the line list will be synthesised.

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

    if 0 >= opacity_contribution:
        raise ValueError("opacity contribution must be a positive float")

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

    if wavelength_region is None:
        wavelength_region = [
            transitions["wavelength"].min() - opacity_contribution,
            transitions["wavelength"].max() + opacity_contribution
        ]

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
        opacity_contribution, in_npoints=pixels, in_modtype=modtype,
        in_debug=debug, #f2pystop=MOOGException(),
        data_path=os.path.dirname(__file__))

    assert wavelengths.size == fluxes.size
    assert wavelengths.size == int(pixels)

    return (wavelengths, fluxes)


def atomic_abundances(transitions, photosphere_information, microturbulence,
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
    safe_mode = kwargs.pop("safe_mode", False)
    modtype, photosphere_arr, metallicity = _format_photosphere(
        photosphere_information, photosphere_kwargs,
        interpolator=kwargs.pop("_interpolator", None))

    # Prepare the transitions table.
    transitions = _format_transitions(transitions)
    
    # Prepare the abundance information
    photospheric_abundances = _format_abundances(photospheric_abundances)

    # Calculate abundances.
    if not safe_mode: 
        code, output = moog.abundances(metallicity, microturbulence,
            photosphere_arr, photospheric_abundances, transitions,
            in_modtype=modtype, in_debug=debug, f2pystop=MOOGException(),
            data_path=os.path.dirname(__file__))
        return output

    else:
        raise NotImplementedError
        print("IN SAFE MODE")
        p = multiprocessing.Process(target=moog.abundances, args=(metallicity,
            microturbulence, photosphere_arr, photospheric_abundances,
            transitions), kwargs={"in_modtype": modtype, "in_debug": debug})
        p.start()
        p.join(10)
        if p.is_alive():
            logger.warn("Terminating MOOG process that has timed out")
            p.terminate()
            p.join()

            return np.ones(len(transitions)) * np.nan
        
        # So. Slow. But now that we know this will work...
        code, output = moog.abundances(metallicity, microturbulence,
            photosphere_arr, photospheric_abundances, transitions,
            in_modtype=modtype, in_debug=debug)
        return output


        


    

