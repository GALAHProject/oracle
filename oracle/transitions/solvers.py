#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" Solvers to fit atomic absorption transitions. """

from __future__ import absolute_import, print_function

__author__ = "Andy Casey <arc@ast.cam.ac.uk>"

import logging
import numpy as np

from scipy import optimize as op
from scipy.special import wofz


logger = logging.getLogger("oracle")


class BaseFitter(object):

    _initialised = False

    def __init__(self, *args, **kwargs):
        self._initialised = True



class gaussian_profile(object):

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
    def integrate(x, mu, sigma, amplitude, **kwargs):
        return amplitude * sigma * np.sqrt(2*np.pi)


class voigt_profile(object):

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
        v = amplitude * np.real(wofz(z)) 
    
        return continuum - v

    @staticmethod
    def integrate(x, mu, sigma, amplitude, shape, continuum=1.0):

        raise NotImplementedError




def mask_data(x, y, mask_regions, mask_non_finites=False):
    """
    Return a mask array for the data (x, y).
    """

    assert x.size == y.size
    if mask_non_finites:
        mask = np.isfinite(y)
    else:
        mask = np.ones(x.size, dtype=bool)

    for lower, upper in mask_regions:
        if None not in (upper, lower):
            assert upper > lower

        mask *= upper > x if upper is not None else 1
        mask *= x > lower if lower is not None else 1

    return mask


class NormalisedSpectrumFitter(BaseFitter):
    """
    Fit absorption profiles to normalised, rest-frame spectra. No accounting for
    continuum or redshift.
    """

    def __init__(self, global_mask=None, profile="gaussian", initial_fwhm=0.05, 
        central_weighting=True, wavelength_tolerance=0.05):
        """
        Fit profiles to absorption lines.
        """

        if profile.lower() not in ("gaussian", "voigt"):
            raise ValueError("profile must be either gaussian or a voigt")

        if 0 >= initial_fwhm:
            raise ValueError("initial FWHM value must be positive")

        if 0 > wavelength_tolerance:
            raise ValueError("wavelength tolerance must be zero or positive")

        self._default_kwargs = {
            "global_mask": global_mask if global_mask is not None else [],
            "profile": profile.lower(),
            "initial_fwhm": initial_fwhm,
            "central_weighting": bool(central_weighting),
            "wavelength_tolerance": wavelength_tolerance,
            "maximum_wavelength_window": 1.0
        }

        return None


    def fit_transition(self, spectrum, wavelength, mask=None, **kwargs):

        # Get the keywords.
        kwds = self._default_kwargs.copy()
        kwds.update(**kwargs)

        if mask is not None:
            # Join the masks so we only have to deal with one.
            kwds["global_mask"].extend(mask)

        if not (spectrum.disp[-1] > wavelength > spectrum.disp[0]):
            raise ValueError("wavelength is outside of spectral range")

        # Get the index of the profile trough.
        if kwds["wavelength_tolerance"] > 0:
            # Look for the lowest flux point between wavelength +/- tolerance.
            idxs = spectrum.disp.searchsorted([
                wavelength - kwds["wavelength_tolerance"],
                wavelength + kwds["wavelength_tolerance"]
            ])
            index = np.nanargmin(spectrum.flux.__getslice__(*idxs)) + idxs[0]
        else:
            index = spectrum.disp.searchsorted(wavelength)

        # What wavelength range will actually be fit?
        indices = spectrum.disp.searchsorted([
            wavelength - kwds["maximum_wavelength_window"],
            wavelength + kwds["maximum_wavelength_window"]
        ])

        disp = spectrum.disp.__getslice__(*indices)
        flux = spectrum.flux.__getslice__(*indices)
        variance = spectrum.variance.__getslice__(*indices)

        # Apply any masks.
        ma = mask_data(disp, flux / variance, kwds["global_mask"],
            mask_non_finites=True)

        if not np.any(ma):
            # No finite data around the line.
            raise ValueError("no finite data within {0:.2f} Angstroms of the "\
                "line at {1:.2f} Angstroms".format(kwds["wavelength_tolerance"],
                    wavelength))


        # Get the profile and initialise the point.
        p_init = np.array([
            spectrum.disp[index],           # mu
            2.355 * kwds["initial_fwhm"],   # sigma
            1.0 - spectrum.flux[index]      # amplitude
        ])
        if kwds["profile"] == "gaussian":
            f = gaussian_profile() 
            
        elif kwds["profile"] == "voigt":
            f = voigt_profile()
            p_init = np.append(p_init, [0.01]) # shape
        

        try:
            # Note: Ensure we are only using finite values for the fit.
            p_opt, p_cov = op.curve_fit(f, disp[ma], flux[ma], p_init,
                sigma=np.sqrt(variance[ma]), absolute_sigma=True)

        except:
            logger.info("Exception occurred during line fitting procedure.")
            raise

        # TODO do something about if the wavelength has moved too much

        fig, ax = plt.subplots()
        ax.plot(disp[ma], flux[ma], c='k')
        ax.plot(disp, f(disp, *p_init), c='g')
        ax.plot(disp, f(disp, *p_opt), c='r')

        # Fit as voigt now
        f = voigt_profile()
        p_init = np.append(p_opt.copy(), [0.0001])

        ax.plot(disp, f(disp, *p_init), "k", lw=2)

        p_opt2, p_cov = op.curve_fit(f, disp[ma], flux[ma], p_init,
            sigma=np.sqrt(variance[ma]), absolute_sigma=True)

        ax.plot(disp, f(disp, *p_opt2), 'r:')

        #f.integrate(disp, *p_opt)

        raise a


        # Fit the profile by least squares.


        # Do a secondary fit by excluding outlier points.

        # Integrate the profile.



if __name__ == "__main__":

    import oracle

    spec = oracle.specutils.Spectrum1D.load_GALAH(
        "/Users/arc/research/galah/data/iDR1/data/benchmark/DeltaEri_1.fits", normalised=True)


    f = NormalisedSpectrumFitter(profile="gaussian")
    f.fit_transition(spec, 4779.29)


