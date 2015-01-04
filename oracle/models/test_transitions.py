
import os

from oracle.specutils import Spectrum1D
import transitions
#from oracle.models import transitions

def atomic_transition():

    data = Spectrum1D.load("Spectra/ccd_3/sun_red.0002.fits")
    atomic_data = {
        "wavelength": 6609.110,
        "species": 26.0,
        "excitation_potential": 2.557,
        "loggf": -2.692,
        "blending_transitions": [
            [6607.652, 24.0, 5.788, -3.784],
            [6607.675, 24.0, 3.365, -4.181],
            [6607.818, 23.0, 1.349, -1.960],
            [6608.025, 26.0, 2.277, -4.030],
            [6608.198, 57.0, 1.234, -0.780],
            [6608.434, 23.0, 4.305, -2.866],
            [6608.473, 21.0, 3.875, -2.298],
            [6608.626, 24.1, 6.482, -3.089],
            [6608.793, 26.0, 5.326, -4.209],
            [6608.936, 24.0, 4.398, -1.048],
            [6608.952, 26.0, 4.188, -2.781],
            [6609.032, 74.0, 1.915, -2.530],
            [6609.234, 26.1, 7.123, -3.022],
            [6609.252, 20.0, 5.736, -5.852],
            [6609.255, 25.1, 6.860, -2.053],
            [6609.678, 26.0, 0.989, -5.002],
            [6609.985, 20.0, 5.737, -4.014],
            [6610.035, 64.1, 1.657, -1.786],
            [6610.073, 20.0, 5.877, -5.517],
            [6610.392, 26.0, 5.503, -2.864],
            [6610.453, 25.0, 5.920, -1.386],
            [6610.475, 21.1, 7.466, 0.325]
        ]
    }

    line = transitions.AtomicTransition(**atomic_data)
    print(line)

    initial_theta = {
        "effective_temperature": 5500,
        "surface_gravity": 4.0,
        "metallicity": 0.,
        "microturbulence": 1.4
    }
    line.fit_profile(data, initial_theta=initial_theta, continuum_order=1)

    raise a

# Don't run this test on Travis.
if "TRAVIS_BUILD_ID" not in os.environ:
    atomic_transition()