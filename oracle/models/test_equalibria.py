
import os
import oracle




def equalibria():
    data = [
        oracle.specutils.Spectrum1D.load("Spectra/ccd_1/sun_blu.0002.fits"),
        oracle.specutils.Spectrum1D.load("Spectra/ccd_2/sun_grn.0002.fits"),
        oracle.specutils.Spectrum1D.load("Spectra/ccd_3/sun_red.0002.fits"),
        oracle.specutils.Spectrum1D.load("Spectra/ccd_4/sun_ir.0002.fits"),
    ]

    model = oracle.models.EqualibriaModel("hermes-atomic-linelist.yaml")

    theta, r_chi_sq, expected_dispersion, expected_flux = model.initial_theta(
        data, full_output=True)

    stellar_parameters = model.estimate_stellar_parameters(data,
        initial_theta=theta)

    raise a

# Don't run this test on Travis.
if "TRAVIS_BUILD_ID" not in os.environ:
    equalibria()