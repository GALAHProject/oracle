#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" A synthesis solver to fit a blended atomic absorption transition. """

from __future__ import absolute_import, print_function

__author__ = "Andy Casey <arc@ast.cam.ac.uk>"

__all__ = ["BaseSynthesisFitter", "SynthesisFitter"]

import logging
import numpy as np

from scipy import (ndimage, optimize as op)
from oracle.transitions import AtomicTransition

try:
    from .base import BaseFitter

except ValueError:
    from base import BaseFitter
    print("TODO NAUGHTY IMPORTY")

logger = logging.getLogger("oracle")


class BaseSynthesisFitter(BaseFitter):
    pass


class SynthesisFitter(BaseSynthesisFitter):
    """
    Fit synthetic spectra to a single, or multiple portions of observed spectra.
    """

    def __init__(self, global_continuum=True, global_resolution=True, **kwargs):
        self._default_kwargs = {
            "mask": [],
            "v_rad_tolerance": 0,
            "continuum_degree": 0
        }
        self._default_kwargs.update(**kwargs)
        self.global_continuum = global_continuum
        self.global_resolution = global_resolution


    def fit(self, spectrum, transitions, synthesiser_factory, **kwargs):
        """
        Fit multiple transitions in a given spectrum. If global_continuum and
        global_resolution are False, this is equivalent to just fitting each
        line individually.

        :param data:
            The observed spectrum.

        :type data:
            :class:`~oracle.specutils.Spectrum1D`

        :param transitions:
            The atomic transition(s) to measure.

        :type transitions:
            list of :class:`~oracle.transitions.AtomicTransition` objects

        :param synthesiser_factory:
            A synthesising factory that takes an (lower, upper) wavelength pair
            and returns a function that accepts one argument (abundance) and
            returns synthetic spectra over the wavelength range (lower, upper).

        :type synthesiser_factory:
            callable
        """

        # Single line only:
        if isinstance(transitions, AtomicTransition):
            return self._fit(spectrum, transitions, synthesiser_factory,
                **kwargs)

        # Multiple lines, but no global behaviour:
        if not self.global_resolution and not self.global_continuum:
            return [self._fit(spectrum, each, synthesiser_factory, **kwargs) \
                for each in transitions]

        # Multiple lines, global behaviour required:
        full_output = kwargs.pop("full_output", False)


        # Possibilities:
        # global continuum (N parameters)
        # global resolution (1 parameter)
        # RV variations on a line-per-line basis
        # abundances for each little region.



        raise NotImplementedError
        return self._fit(*args, **kwargs)



    def _fit(self, spectrum, transition, synthesiser_factory, **kwargs):
        """
        Fit a portion of spectrum by solving for the abundance, continuum, 
        spectral resolution, and small shifts in radial velocity. This class is
        well-suited to fitting an atomic absorption feature that is blended.

        :param data:
            The observed spectrum.

        :type data:
            :class:`~oracle.specutils.Spectrum1D`

        :param transition:
            The atomic transition to measure.

        :type transition:
            :class:`~oracle.transitions.AtomicTransition`

        :param synthesiser_factory:
            A synthesising factory that takes an (lower, upper) wavelength pair
            and returns a function that accepts one argument (abundance) and
            returns synthetic spectra over the wavelength range (lower, upper).

        :type synthesiser_factory:
            callable
        """

        if not isinstance(transition, AtomicTransition):
            raise TypeError("transition is expected to be an oracle.transitions"
                ".AtomicTransition")

        full_output = kwargs.pop("full_output", False)

        # Get the keywords and update the behaviour keywords.
        kwds = self._default_kwargs.copy()
        kwds["mask"].extend([] if transition.mask is None else transition.mask)
        kwds["continuum_degree"] = transition.continuum_degree
        kwds["v_rad_tolerance"] = transition.v_rad_tolerance
        
        # Set up the synthesiser_factory for this wavelength range.
        synthesiser = synthesiser_factory(transition.fitting_region)
        indices = spectrum.disp.searchsorted(transition.fitting_region)

        disp = spectrum.disp.__getslice__(*indices)
        flux = spectrum.flux.__getslice__(*indices)
        variance = spectrum.variance.__getslice__(*indices)

        # Create a mask.
        mask = self.mask_data(disp, flux / variance, kwds["mask"],
            mask_non_finites=True)

        # Establish initial guess of parameters.
        # What initial parameters could we have?
        # - abundance (always)
        # - spectral resolution (always)
        # - radial velocity (sometimes)
        # - continuum coefficients (sometimes)

        p_init = kwargs.get("p_init", None)
        if p_init is None:
            p_init = np.array([
                0,      # Abundance enhancement releative to photosphere
                20000,  # Spectral resolving power (typical for abundance work)
            ])

            # Any radial velocity?
            if kwds["v_rad_tolerance"] > 0:
                p_init = np.append(p_init, [0])

            # Any continuum?
            continuum_degree = kwds["continuum_degree"]
            if continuum_degree > 0:
                coefficients = np.zeros(continuum_degree)
                coefficients[-1] = np.nanmedian(flux)
                p_init = np.append(p_init, coefficients)

        # Create a model for the data.
        def model(disp, abundance, R, *args):

            _error = np.inf * np.ones_like(disp)

            # If *args is present, infer whether it contains radial velocity,
            # continuum coefficients, or both.
            v_rad, continuum_coefficients = 0, [1]
            if len(args) > 0:
                if kwds["v_rad_tolerance"]:
                    v_rad = args[0]
                    continuum_coefficients = args[1:] if len(args) > 1 else [1]
                else:
                    continuum_coefficients = [] + args

            if (kwds["v_rad_tolerance"] > 0 \
            and abs(v_rad) > kwds["v_rad_tolerance"]) \
            or 0 >= R:
                return _error

            # Synthesise spectra for this abundance.
            try:
                model_disp, model_intensities = synthesiser(abundance)

            except ValueError:
                logger.exception("Exception raised when synthesising region.")
                return _error

            # Apply transformations (radial velocity, continuum, smoothing)
            model_disp *= 1. + v_rad/299792.458
            sigma = model_disp.mean()/(2.35482 * R * np.diff(model_disp).mean())
            model_fluxes = np.polyval(continuum_coefficients, model_disp) \
                * ndimage.gaussian_filter(model_intensities, sigma, cval=np.nan)

            # Interpolate model fluxes onto the observed pixels.
            model_pixels = np.interp(disp, model_disp, model_fluxes,
                left=np.nan, right=np.nan)
            return model_pixels

        try:
            p_opt, p_cov = op.curve_fit(model, disp[~mask], flux[~mask], p_init,
                sigma=np.sqrt(variance[~mask]), absolute_sigma=True, epsfcn=0.0,
                ftol=1e-10, gtol=1e-10)

        except:
            logger.exception("Exception occurred during synthesis fitting:")
            raise

        # Calculate the chi-squared value.
        difference = model(disp[~mask], *p_opt) - flux[~mask]
        chi_sq = sum(difference**2/variance[~mask])
        dof = sum(~mask) - len(p_opt) - 1

        if full_output:
            return (p_opt, chi_sq, dof, kwds,
                [
                    disp[~mask],
                    model(disp[~mask], *p_opt),
                    model(disp[~mask], *p_init)
                ],
                p_cov)

        return (p_opt, chi_sq, dof, kwds)



if __name__ == "__main__":

    # A test case for the synthesis fitter.

    import oracle
    data = oracle.specutils.Spectrum1D.load_GALAH(
        "/Users/arc/research/galah/data/iDR1/data/benchmark/18Sco_3.fits", normalised=True, rest=True)

    # Create a synthesiser that includes the line list and model photosphere
    # information for the fitter.
    interpolator = oracle.photospheres.interpolator(kind="MARCS")
    photosphere = interpolator(5810, 4.45, 0)

    from astropy.table import Table
    t = Table(rows=[{
        "wavelength": 6592.9124, 
        "species": 26.0,
        "excitation_potential": 2.727,
        "loggf": -1.473,
    }])

    import cPickle as pickle
    with open("lines.pkl", "rb") as fp:
        rows = pickle.load(fp)

    for row in rows:
        row["VDW_DAMP"] = row["C6"]
    t = Table(rows=rows)
    
    synthesiser_factory = lambda wavelength_range: \
        lambda abundance: oracle.synthesis.moog.synthesise(t,
            photosphere, wavelength_range, microturbulence=1.07,
            photospheric_abundances=[26, abundance], damping=4)

    # 0 = 840699499.20031548
    # 1 = 111950631.80524081
    # 2 = 126757416.63107705
    # 3 = 1231810711.2760828
    # 4 = 116284148.07045062


    foo = SynthesisFitter()

    transition = oracle.transitions.AtomicTransition(wavelength=6592.9124, 
        species=26.0, e_low=2.727, log_gf=-1.473)
#    transition.mask = [[6593.44, 6594.27]]

    p_opt, chi_sq, dof, kwds, (disp, opt_mod_flux, init_mod_flux), cov \
        = foo.fit(data, transition, synthesiser_factory,
            full_output=True)


    import matplotlib.pyplot as plt
    fig, ax = plt.subplots()
    ax.plot(disp, init_mod_flux, c='b', zorder=10)
    ax.plot(disp, opt_mod_flux, c='r', zorder=10)
    lims = ax.get_xlim()
    
    ax.errorbar(data.disp, data.flux, yerr=data.variance**0.5, fmt=None, ecolor="k")

    ax.plot(data.disp, data.flux, c='k')
    ax.set_xlim(lims)

    raise a


