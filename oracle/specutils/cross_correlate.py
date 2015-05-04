# coding: utf-8

""" Cross-correlating spectra. """

from __future__ import division, print_function

__author__ = "Andy Casey <arc@ast.cam.ac.uk>"

__all__ = ["cross_correlate", "cross_correlate_grid"]

import numpy as np
import multiprocessing
import astropy.units as u
from astropy import modeling, constants


def cross_correlate(observed, template, wavelength_range=None):
    """
    Return a redshift by cross correlation of a template and observed spectra.

    :param observed:
        The observed spectrum.

    :type observed:
        :class:`Spectrum1D`

    :param template:
        The template spectrum, expected to be at rest-frame.

    :type template:
        :class:`Spectrum1D`

    :param wavelength_range: [optional]
        The (observed) start and ending wavelengths to use for the cross correlation.

    :type wavelength_range:
        tuple

    :returns:
        The relative velocity and associated uncertainty in km/s.

    :rtype:
        tuple
    """

    # Put the spectra on the same dispersion mapping
    if wavelength_range is not None:
        if not isinstance(wavelength_range, (list, tuple, np.array)) \
        or len(wavelength_range) != 2:
            raise TypeError("wavelength range must either be None or a two-length"\
                " tuple-like object with the start and end wavelength ranges")

        indices = observed.disp.searchsorted(wavelength_range)
        dispersion = observed.disp[indices[0]:indices[1] + 1]
        observed_flux = observed.flux[indices[0]:indices[1] + 1]

    else:
        dispersion = observed.disp
        observed_flux = observed.flux

    template_flux = np.interp(dispersion, template.disp, template.flux,
                              left=1, right=1)

    # Be forgiving, although we shouldn't have to be.
    N = np.min(map(len, [dispersion, observed_flux, template_flux]))

    # Ensure an even number of points
    if N % 2 > 0:
        N -= 1

    dispersion = dispersion[:N]
    observed_flux = observed_flux[:N]
    template_flux = template_flux[:N]

    assert len(dispersion) == len(observed_flux)
    assert len(observed_flux) == len(template_flux)

    # Set up z array
    m = len(dispersion) / 2
    z_array = dispersion/dispersion[N/2] - 1.0

    # Apodize edges
    edge_buffer = 0.1 * (dispersion[-1] - dispersion[0])
    low_w_indices = np.nonzero(dispersion < dispersion[0] + edge_buffer)[0]
    high_w_indices = np.nonzero(dispersion > dispersion[-1] - edge_buffer)[0]

    apod_curve = np.ones(N, dtype='d')
    apod_curve[low_w_indices] = (1.0 + np.cos(np.pi*(1.0 - \
        (dispersion[low_w_indices] - dispersion[0])/edge_buffer)))/2.
    apod_curve[high_w_indices] = (1.0 + np.cos(np.pi*(1.0 - \
        (dispersion[-1] - dispersion[high_w_indices])/edge_buffer)))/2.

    apod_observed_flux = observed_flux * apod_curve
    apod_template_flux = template_flux * apod_curve

    fft_observed_flux = np.fft.fft(apod_observed_flux)
    fft_template_flux = np.fft.fft(apod_template_flux)
    template_flux_corr = (fft_observed_flux * fft_template_flux.conjugate())   \
        / np.sqrt(np.inner(apod_observed_flux, apod_observed_flux)             \
        * np.inner(apod_template_flux, apod_template_flux))

    correlation = np.fft.ifft(template_flux_corr).real

    # Reflect about zero
    ccf = np.zeros(N)
    ccf[:N/2] = correlation[N/2:]
    ccf[N/2:] = correlation[:N/2]

    # Get height and redshift of best peak
    # TODO: Fit a Gaussian profile here instead, and get the z_err from the FWHM
    #       of the profile

    # Scale the CCF
    h = ccf.max()
    ccf -= ccf.min()
    ccf *= (h/ccf.max())

    initial_stddev = np.ptp(z_array[np.where(ccf >= 0.5*h)])

    model = modeling.models.Gaussian1D(
        mean=z_array[np.argmax(ccf)], amplitude=1,
        stddev=0.05)

    fitter = modeling.fitting.LevMarLSQFitter()

    profile = fitter(model, z_array, ccf)

    import matplotlib.pyplot as plt
    plt.plot(z_array, ccf, 'k')
    plt.plot(z_array, profile(z_array), 'r')

    raise a



    #c = 299792.458 # km/s
    #z_best = z_array[ccf.argmax()]
    #z_err = /2.35482)**2

    return (z_best * c, z_err * c)


def _ccf(z_array, apod_template_flux, template_flux_corr, N, index=None,
         method="fast", **kwargs):

    denominator = np.sqrt(np.inner(apod_template_flux, apod_template_flux))
    flux_correlation = template_flux_corr / denominator
    correlation = np.fft.ifft(flux_correlation).real

    # Reflect about zero
    ccf = np.zeros(N)
    ccf[:N/2] = correlation[N/2:]
    ccf[N/2:] = correlation[:N/2]

    # Get height and redshift of best peak
    h = ccf.max()

    # Scale the CCF
    ccf -= ccf.min()
    ccf *= (h/ccf.max())

    if method != "fast":
        mean = ccf.argmax()

        model = modeling.models.Gaussian1D(mean=mean, amplitude=h,
                                           stddev=1)
        model += modeling.models.Const1D()

        fitter = modeling.fitting.LevMarLSQFitter()

        fit_sigma = kwargs.pop("__fit_sigma", 1.5)
        x = np.arange(N)
        # TODO revise this range?
        use = (mean + 5 >= x) * (x >= mean - 5)
        model.bounds["stddev_0"] = [0, 1.5]
        profile = fitter(model, x[use], ccf[use])

        mean, stddev, amplitude = profile.mean_0.value, \
            profile.stddev_0.value * np.diff(z_array)[0], \
            profile.amplitude_0.value + profile.amplitude_1.value

    else:

        stddev = np.ptp(z_array[np.where(ccf >= 0.5 * h)]/2.355)
        mean, amplitude = ccf.argmax(), h

    mean = np.interp(mean, np.arange(N), z_array)

    if index is not None:
        return (index, mean, stddev, amplitude)
    return (mean, stddev, amplitude)


def cross_correlate_grid(template_dispersion, template_fluxes, observed_flux,
                         continuum_degree=4, apodize=0.10, remeasure_best=True,
                         threads=1):

    if template_dispersion.shape[0] != template_fluxes.shape[1]:
        raise ValueError("template dispersion must have size (N_pixels,) and "\
            "template fluxes must have size (N_models, N_pixels)")
    try:
        continuum_degree = int(continuum_degree + 1)
    except (TypeError, ValueError):
        raise TypeError("continuum order must be an integer-like object")

    assert 1 > apodize >= 0, "Apodisation fraction must be between 0 and 1"

    N = template_dispersion.size
    N = N - 1 if N % 2 > 0 else N
    N_models = template_fluxes.shape[0]

    dispersion = template_dispersion[:N]
    template_flux = template_fluxes[:, :N]

    observed_flux = observed_flux.copy()[:N]
    non_finite = ~np.isfinite(observed_flux)
    if non_finite.sum() > 0:
        observed_flux[non_finite] = np.interp(dispersion[non_finite],
            dispersion[~non_finite], observed_flux[~non_finite])

    # Normalise
    if continuum_degree >= 1:
        coeffs = np.polyfit(dispersion, observed_flux, continuum_degree)
        observed_flux /= np.polyval(coeffs, dispersion)

    # Scale the flux level to that the template intensities
    observed_flux = (observed_flux * template_flux.ptp()) + template_flux.min()

    # Apodize edges
    edge_buffer = apodize * (dispersion[-1] - dispersion[0])
    low_w_indices = np.nonzero(dispersion < dispersion[0] + edge_buffer)[0]
    high_w_indices = np.nonzero(dispersion > dispersion[-1] - edge_buffer)[0]

    apod_curve = np.ones(N, dtype='d')
    apod_curve[low_w_indices] = (1.0 + np.cos(np.pi*(
        1.0 - (dispersion[low_w_indices] - dispersion[0])/edge_buffer)))/2.
    apod_curve[high_w_indices] = (1.0 + np.cos(np.pi*(
        1.0 - (dispersion[-1] - dispersion[high_w_indices])/edge_buffer)))/2.

    apod_observed_flux = observed_flux * apod_curve
    apod_template_flux = template_flux * apod_curve

    fft_observed_flux = np.fft.fft(apod_observed_flux)
    fft_template_flux = np.fft.fft(apod_template_flux)
    template_flux_corr = (fft_observed_flux * fft_template_flux.conjugate())
    template_flux_corr /= np.sqrt(np.inner(apod_observed_flux, apod_observed_flux))

    z_array = np.array(dispersion.copy())/dispersion[N/2] - 1.0

    z = np.ones(N_models) * np.nan
    z_err = np.ones(N_models) * np.nan
    R = np.ones(N_models) * np.nan

    if threads > 1:
        processes = []
        pool = multiprocessing.Pool(threads)
        for i in xrange(N_models):
            processes.append(pool.apply_async(_ccf, args=(z_array,
                apod_template_flux[i, :], template_flux_corr[i, :], N, i)))

        for process in processes:
            result = process.get()
            index = result[0]
            z[index], z_err[index], R[index] = result[1:]

        pool.close()
        pool.join()

    else:
        # Thread this!
        for i in xrange(N_models):
            z[i], z_err[i], R[i] = _ccf(z_array, apod_template_flux[i, :],
                template_flux_corr[i, :], N)


    c = constants.c.to("km/s").value

    # Should we precisely re-measure the best point?
    if remeasure_best:
        best = R.argmax()

        b_z, b_z_err, b_R = _ccf(z_array, apod_observed_flux,
            template_flux_corr[best, :], N, None, method="slow")

        z[best] = b_z
        z_err[best] = b_z_err

    return (z * c, z_err * c, R)
