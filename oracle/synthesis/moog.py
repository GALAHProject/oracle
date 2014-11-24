

import numpy as np

import _moog
import oracle.atmospheres



def synthesize(teff, logg):

    photospheric_structure = oracle.atmospheres.interpolate(teff, logg, metallicity,
        **atmosphere_kwargs)

    return None


def abundances(teff, logg, metallicity, microturbulence, transitions,
    photospheric_abundances=None, atmosphere_kwargs=None):

    # Transitions should be an array of minimum size, or a record array to
    # indicate what is what

    # Interpolate the photospheric structure
    photospheric_structure = oracle.atmospheres.interpolate(teff, logg,
        metallicity, **atmosphere_kwargs)


    if photospheric_abundances is None:
        photospheric_abundances = np.asfortranarray(np.array([])).reshape(-1, 2)
    else:
        photospheric_abundances = np.asfortranarray(
            np.atleast_2d(photospheric_abundances))

    code, output = _moog.abundances(teff, logg, metallicity, microturbulence,
        photospheric_structure, photospheric_abundances, transitions)

    return output