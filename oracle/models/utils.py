#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" Utility functions that might be refactored later. """

from __future__ import absolute_import, print_function

__author__ = "Andy Casey <arc@ast.cam.ac.uk>"

import logging
import numpy as np
from scipy import stats

from oracle.atmospheres import solar_abundance

logger = logging.getLogger(__name__)

def equalibrium_state(data, sigma_clip=5, metallicity=None, full_output=False):
    """
    Calculate the current equalibrium state.
    """

    # Use the abundances in self.atomic_transitions
    finite = np.isfinite(data["abundance"])
    if finite.sum() == 0:
        raise ValueError("no abundances calculated")

    # Slice out just the lines with abundances.
    data = data[finite]

    outliers_removed, outliers = False, np.zeros(len(data), dtype=bool)
    neutral = ((data["species"] % 1) == 0)
    ionised = ((data["species"] % 1) > 0)
    if metallicity is None:
        metallicity = data.meta["abundances_given_stellar_parameters"][2]

    while True:
        neutral *= ~outliers
        ionised *= ~outliers
        combined = neutral + ionised

        # Calculate the slope with excitation potential and abundance.
        excitation_regression = stats.linregress(
            x=data["excitation_potential"][combined],
            y=data["abundance"][combined])

        # Calculate the ionisation state.
        ionisation_state = np.median(data["abundance"][neutral]) \
            - np.median(data["abundance"][ionised])

        # Calculate the abundance state.
        abundance_state = np.median(data["abundance"][combined] - \
            (solar_abundance(data["species"][combined]) + metallicity))

        # Slope with reduced equivalent width and line abundance.
        rew = np.log(data["equivalent_width"] / data["wavelength"])
        rew_regression = stats.linregress(
            x=rew[combined],
            y=data["abundance"][combined])

        # Calculate the state
        state = np.array([
            excitation_regression[0],
            ionisation_state,
            abundance_state,
            rew_regression[0]
        ])

        if outliers_removed:
            final_state = (excitation_regression, ionisation_state,
                abundance_state, rew_regression)
            break

        else:
            initial_state = (excitation_regression, ionisation_state,
                abundance_state, rew_regression)

            if sigma_clip is None or 0 >= sigma_clip:
                # Don't remove any outliers
                outliers_removed, final_state = False, initial_state
                break

            # We don't want to remove stars that are simply on the edge of
            # distributions because if our initial guess was very wrong, we
            # could just be removing the weakest or strongest lines.

            # Instead we will remove lines that are discrepant from the line
            # fits.

            # Outliers in excitation potential vs abundance:
            x = data["excitation_potential"][combined]
            y = data["abundance"][combined]
            line = excitation_regression[0] * x + excitation_regression[1]

            differences = np.abs(line - y)
            excitation_sigma = differences/np.std(differences)

            # Outliers in reduced equivalent width vs abundance:
            x = rew[combined]
            y = data["abundance"][combined]
            line = rew_regression[0] * x + rew_regression[1]

            differences = np.abs(line - y)
            line_strength_sigma = differences/np.std(differences)

            # Update the finite mask to remove outliers
            logger.debug("Largest deviant was at {0:.2f} sigma".format(
                max([line_strength_sigma.max(), excitation_sigma.max()])))


            outliers_to_remove = (excitation_sigma > sigma_clip) \
                + (line_strength_sigma > sigma_clip)

            # Check that we aren't throwing out the last of a species.
            if len(set(data["species"][combined][~outliers_to_remove])) < \
                len(set(data["species"][combined])):
                
                species_being_removed = list(set(data["species"][combined])\
                    .difference(data["species"][combined][~outliers_to_remove]))

                logger.warn("Doing outlier rejection would remove the last {0} "
                    "line. Instead we will just remove the most discrepant"
                    "measurement from that species, then hope.".format(
                        species_being_removed[0]))

                _ = (data["species"][combined] == species_being_removed[0])
                __ = np.argsort((excitation_sigma + line_strength_sigma)[np.where(_)[0]])[0]
                outliers_to_remove[np.where(_)[0][__]] = False

            outliers[combined] = outliers_to_remove
            outliers_removed = True
            continue # to re-fit the lines

    if full_output:
        info = {
            "initial_state": initial_state,
            "final_state": final_state,
            "coefficients": [
                [excitation_regression[0], excitation_regression[1]],
                [rew_regression[0], rew_regression[1]]
            ],
            "outliers": np.where(finite)[0][outliers],
            "~outliers": np.where(finite)[0][~outliers]
        }
        return (state, info)
    return state


def jacobian_original(stellar_parameters, *args, **kwargs):

    teff, vt, logg, feh = stellar_parameters

    # This is the black magic.
    full_jacobian = np.array([
        [ 5.4393e-08*teff - 4.8623e-04, -7.2560e-02*vt + 1.2853e-01,  1.6258e-02*logg - 8.2654e-02,  1.0897e-02*feh - 2.3837e-02],
        [ 4.2613e-08*teff - 4.2039e-04, -4.3985e-01*vt + 8.0592e-02, -5.7948e-02*logg - 1.2402e-01, -1.1533e-01*feh - 9.2341e-02],
        [-3.2710e-08*teff + 2.8178e-04,  3.8185e-03*vt - 1.6601e-02, -1.2006e-02*logg - 3.5816e-03, -2.8592e-05*feh + 1.4257e-03],
        [-1.7822e-08*teff + 1.8250e-04,  3.5564e-02*vt - 1.1024e-01, -1.2114e-02*logg + 4.1779e-02, -1.8847e-02*feh - 1.0949e-01]
    ])
    return full_jacobian.T


def jacobian(stellar_parameters, *args, **kwargs):
    """
    Calculate the approximate Jacobian matrix, given some stellar parameters.

    :param stellar_parameters:
        The stellar parameters to calculate the Jacobian at. These are expected
        to be the effective temperature, the surface gravity, the overall
        metallicity, and the microturbulence.

    :type stellar_parameters:
        :class:`numpy.array`

    :returns:
        The Jacobian approximation.
    """

    #logger.debug("Updating Jacobian at Teff = {0:.0f} K, logg = {1:.2f}, [M/H] "
    #    "= {2:.2f}, vt = {3:.2f}".format(*stellar_parameters))

    # (effective_temperature, surface_gravity, metallicity, microturbulence)
    m = np.array([
        [ 5.4393e-08,  1.6258e-02,  1.0897e-02, -7.2560e-02],
        [ 4.2613e-08, -5.7948e-02, -1.1533e-01, -4.3985e-01],
        [-3.2710e-08, -1.2006e-02, -2.8592e-05,  3.8185e-03],
        [-1.7822e-08, -1.2114e-02, -1.8847e-02,  3.5564e-02]
    ])
    b = np.array([
        [-4.8623e-04, -8.2654e-02, -2.3837e-02, +1.2853e-01],
        [-4.2039e-04, -1.2402e-01, -9.2341e-02, +8.0592e-02],
        [+2.8178e-04, -3.5816e-03, +1.4257e-03, -1.6601e-02],
        [+1.8250e-04, +4.1779e-02, -1.0949e-01, -1.1024e-01]
    ])

    n = stellar_parameters.size
    x = np.repeat(stellar_parameters, n).reshape(n, n).T
    jacobian = (m * x + b).T
    x1 = stellar_parameters + np.dot(stellar_parameters,jacobian)
    logger.debug("Jacobian {0} predicts {1}".format(stellar_parameters, x1))

    return jacobian