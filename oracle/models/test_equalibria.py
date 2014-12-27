
import oracle

data = [oracle.specutils.Spectrum1D.load("Spectra/ccd_1/sun_blu.0002.fits")]

model = oracle.models.EqualibriaModel({
        "model": {
            "redshift": True,
            "instrumental_resolution": True,
            "continuum": 2,
            "atomic_lines": [
                # wavelength, species, chi, loggf, damp1, damp2, synthesise_surrounding, opacity_contribution
                [4800.648, 26.0, 4.120, -1.028, 0, 0, 1.0],
                [4871.318, 26.0, 2.870, -0.362]
            ]
#            "atomic_lines_filename": "test_atomic_lines.txt"
#            "blending_lines_filename": "test_blen"
        },
        "settings": {
            "threads": 4
        }
    })


theta = model.initial_theta(data)
