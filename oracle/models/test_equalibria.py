
import os
import oracle




def equalibria():
    data = [
        oracle.specutils.Spectrum1D.load("Spectra/ccd_1/sun_blu.0002.fits"),
        oracle.specutils.Spectrum1D.load("Spectra/ccd_2/sun_grn.0002.fits"),
        oracle.specutils.Spectrum1D.load("Spectra/ccd_3/sun_red.0002.fits"),
        oracle.specutils.Spectrum1D.load("Spectra/ccd_4/sun_ir.0002.fits"),
    ]

    data = [
        oracle.specutils.Spectrum1D.load("/Users/arc/codes/oracle-with-siu/oracle/tests/data/benchmarks/18Sco/18Sco_narval_blue_noresample.txt"),
        oracle.specutils.Spectrum1D.load("/Users/arc/codes/oracle-with-siu/oracle/tests/data/benchmarks/18Sco/18Sco_narval_green_noresample.txt"),
        oracle.specutils.Spectrum1D.load("/Users/arc/codes/oracle-with-siu/oracle/tests/data/benchmarks/18Sco/18Sco_narval_red_noresample.txt"),
        oracle.specutils.Spectrum1D.load("/Users/arc/codes/oracle-with-siu/oracle/tests/data/benchmarks/18Sco/18Sco_narval_ir_noresample.txt")
    ]

    model = oracle.models.EqualibriaModel("hermes-classical.yaml")

    initial_theta, r_chi_sq, expected_dispersion, expected_flux = model.initial_theta(
        data, full_output=True)

    fig, axes = plt.subplots(4)
    indices = [-1] + list(np.where(np.diff(expected_dispersion) > 10)[0]) + [None]
    for i, ax in enumerate(axes):

        disp = expected_dispersion[indices[i]+1:indices[i+1]]
        flux = expected_flux[indices[i]+1:indices[i+1]]

        ax.plot(disp, flux, "bgrr"[i])
        ax.plot(data[i].disp, data[i].flux, 'k')


    stellar_parameters = model.estimate_stellar_parameters(data,
        initial_theta=initial_theta)

    raise a

# Don't run this test on Travis.
if "TRAVIS_BUILD_ID" not in os.environ:
    equalibria()