#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" A synthesis solver to fit a blended atomic absorption transition. """

from __future__ import absolute_import, print_function

__author__ = "Andy Casey <arc@ast.cam.ac.uk>"

import logging
import numpy as np

from scipy import (ndimage, optimize as op)
try:
    from .base import BaseFitter

except ValueError:
    from base import BaseFitter
    print("TODO NAUGHTY IMPORTY")

logger = logging.getLogger("oracle")

class SynthesisFitter(BaseFitter):
    """
    Fit synthetic spectra to a small portion of spectrum. This class is intended
    for fitting a single blended absorption line.
    """

    def __init__(self, global_mask=None, radial_velocity_tolerance=0,
        continuum_degree=0):

        self._default_kwargs = {
            "global_mask": global_mask if global_mask is not None else [],
            "radial_velocity_tolerance": radial_velocity_tolerance,
            "continuum_degree": max(0, continuum_degree)
        }


    def fit(self, spectrum, synthesiser, mask=None, **kwargs):
        """
        Fit a portion of spectrum by solving for the abundance, continuum, 
        spectral resolution, and small shifts in radial velocity. This class is
        well-suited to solving for atomic absorption features that are blended.

        :param data:
            The observed spectrum.

        :type data:
            :class:`~oracle.specutils.Spectrum1D`

        :param synthesiser:
            A synthesising function that will take an abundance and produce a 
            spectrum over a small wavelength region. Thus, this fitting routine
            is subject to the synthesiser class, which assumes the stellar
            parameters, line list and wavelength region remain constant.

        """

        full_output = kwargs.pop("full_output", False)

        # Get the keywords and update the mask if provided.
        kwds = self._default_kwargs.copy()
        kwds.update(**kwargs)
        if mask is not None:
            kwds["global_mask"].extend(mask)

        # What initial parameters could we have?
        # - abundance (always)
        # - spectral resolution (always)
        # - radial velocity (sometimes)
        # - continuum coefficients (sometimes)
        
        # Get the wavelength range synthesised by doing a test synthesis call.
        try:
            d, _ = synthesiser(0)

        except:
            logger.exception("Exception raised when trying synthesiser")
            raise

        indices = spectrum.disp.searchsorted([d[0], d[-1]])

        disp = spectrum.disp.__getslice__(*indices)
        flux = spectrum.flux.__getslice__(*indices)
        variance = spectrum.variance.__getslice__(*indices)

        # Create a mask.
        mask = self.mask_data(disp, flux / variance, kwds["global_mask"],
            mask_non_finites=True)

        # Establish initial guess of parameters.
        p_init = kwargs.get("p_init", None)
        if p_init is None:
            p_init = np.array([
                0,      # Abundance enhancement releative to photosphere
                20000,  # Spectral resolving power (typical for abundance work)
            ])

            # Any radial velocity?
            if kwds["radial_velocity_tolerance"] > 0:
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
                if kwds["radial_velocity_tolerance"]:
                    v_rad = args[0]
                    continuum_coefficients = args[1:] if len(args) > 1 else [1]
                else:
                    continuum_coefficients = [] + args

            if (kwds["radial_velocity_tolerance"] > 0 \
            and abs(v_rad) > kwds["radial_velocity_tolerance"]) \
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

        # TODO: Do we want to do sigma-clipping and go again? Probably not...

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
                ], p_cov)

        return (p_opt, chi_sq, dof, kwds)



if __name__ == "__main__":

    # A test case for the synthesis fitter.

    import oracle
    data = oracle.specutils.Spectrum1D.load_GALAH(
        "/Users/arc/research/galah/data/iDR1/data/benchmark/18Sco_3.fits", normalised=True, rest=True)

    # Create a synthesiser that includes the line list and model atmosphere
    # information for the fitter.
    interpolator = oracle.atmospheres.interpolator(kind="MARCS")
    photosphere = interpolator(5810, 4.45, 0)

    from astropy.table import Table
    transitions = Table(rows=[{
        "wavelength": 6592.9124, 
        "species": 26.0,
        "excitation_potential": 2.727,
        "loggf": -1.473,
        "C6": 320.26400756
    }])

    import cPickle as pickle
    with open("lines.pkl", "rb") as fp:
        rows = pickle.load(fp)

    #transitions = Table(rows=rows)

    wavelength_range = [6592.6880, 6593.4500]
    wavelength_range = [6591.2880, 6594.6500]
    #wavelength_range = [6592.9124 - 5., 6597.2]
    
    synthesiser = lambda abundance: oracle.synthesis.moog.synthesise(transitions,
        photosphere, wavelength_range, microturbulence=1.07, photospheric_abundances=[26, abundance])


    foo = SynthesisFitter(radial_velocity_tolerance=3, continuum_degree=0)
    p_opt, chi_sq, dof, kwds, (disp, opt_mod_flux, init_mod_flux), cov \
        = foo.fit(data, synthesiser, mask=[[6593.44, 6594.27]],
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


