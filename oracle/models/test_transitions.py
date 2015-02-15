
import os
import yaml

import numpy as np

from oracle.specutils import Spectrum1D
import transitions
#from oracle.models import transitions


import matplotlib.pyplot as plt


def minimum_pixel_sampling(data, wavelength):
    """
    Find the minimum pixel size in the data at the given wavelength.
    """

    pixel_size = np.nan
    channel_indices = []
    for i, spectrum in enumerate(data):
        if spectrum.disp[-1] >= wavelength >= spectrum.disp[0]:
            channel_indices.append(i)

            index = spectrum.disp.searchsorted(wavelength)

            diff = np.diff(spectrum.disp)
            if pixel_size > diff[index] or not np.isfinite(pixel_size):
                pixel_size = diff[index]

    if len(channel_indices) == 0 and not np.isfinite(pixel_size):
        raise ValueError("cannot find wavelength {0:.2f} in any data channel"\
            .format(wavelength))

    return (pixel_size, channel_indices)


def atomic_transition():

    datas = [
        Spectrum1D.load("Spectra/ccd_1/sun_blu.0002.fits"),
        Spectrum1D.load("Spectra/ccd_2/sun_grn.0002.fits"),
        Spectrum1D.load("Spectra/ccd_3/sun_red.0002.fits"),
        Spectrum1D.load("Spectra/ccd_4/sun_ir.0002.fits"),
    ]
    with open("hermes-atomic-linelist.yaml", "r") as fp:
        model = yaml.load(fp)

    atomic_transitions = model["model"]["atomic_transitions"]

    # Do all
    failed = 0
    for atomic_data in atomic_transitions:

        try:
            ps, i = minimum_pixel_sampling(datas, atomic_data["wavelength"])
        except ValueError:
            continue
        data = datas[i[0]]

        line = transitions.AtomicTransition(**atomic_data)
        print(line)

        initial_theta = {
            "effective_temperature": 5500,
            "surface_gravity": 4.0,
            "metallicity": 0.,
            "microturbulence": 1.4
        }
        try:
            result = line.fit_profile(data, initial_theta=initial_theta,
                constrain_parameters={"fwhm": [0, 1], "blending_fwhm": [0, 1],
                    "wavelength": [line.wavelength - 0.05, line.wavelength + 0.05]},
                continuum_order=1, full_output=True)
        except ValueError:
            failed += 1
            continue

        fig, ax = plt.subplots()
        ax.plot(result[3].disp, result[3].flux, 'b')
        ax.plot(data.disp, data.flux, 'k')
        ax.set_xlim(result[3].disp[0], result[3].disp[-1])


    raise a




# Don't run this test on Travis.
if "TRAVIS_BUILD_ID" not in os.environ:
    atomic_transition()