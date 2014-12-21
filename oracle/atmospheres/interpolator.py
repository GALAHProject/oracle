# coding: utf-8

""" Interpolate model atmospheres (MARCS only, at this stage) """

from __future__ import division, absolute_import, print_function

import cPickle as pickle
import numpy as np
import scipy.interpolate
from pkg_resources import resource_stream


class Interpolator(object):
    """
    A MARCS model interpolator class. This interpolation routine finds the
    nearest opposing neighbours in all dimensions, rescales all photospheric
    quantities onto a common opacity scale (tau(5000)) using a spline function,
    then interpolates the photospheric quantities using the scaling constants
    calculated by Thomas Masseron (Masseron, 2006).
    """

    _coefficients = {
        "T": [0.15, 0.3, "1-(teff/4000)**2"],
        "logPe": [0.15, 0.06, "1-(teff/3500)**2.5"],
        "logPg": [-0.4, 0.06, "1-(teff/4100)**4"],
    }

    def __init__(self, pickled_atmospheres=None):
        """
        Loads the pickled atmospheres.
        """

        if pickled_atmospheres is None:
            with resource_stream(__name__, "marcs-2011-standard.pickle") as fp:
                stellar_parameters, photospheres, photosphere_cols, comment \
                    = pickle.load(fp)
        else:
            with open(pickled_atmospheres, "rb") as fp:
                stellar_parameters, photospheres, photosphere_cols, comment \
                    = pickle.load(fp)

        # Look for duplicate stellar parameter rows
        array_view = stellar_parameters.view(float).reshape(
            stellar_parameters.size, -1)
        _ = np.ascontiguousarray(array_view).view(np.dtype((np.void,
            array_view.dtype.itemsize * array_view.shape[1])))
        _, idx = np.unique(_, return_index=True)

        if idx.size != stellar_parameters.size:
            raise ValueError("{} duplicate stellar parameters found".format(
                stellar_parameters.size - idx.size))

        self.stellar_parameters = stellar_parameters
        self.photospheres = photospheres
        self.photosphere_cols = photosphere_cols
        self.comment = comment

        # Create unique copies of stellar parameters for faster access
        names = stellar_parameters.dtype.names
        self._stellar_parameters = dict(zip(names,
            [np.unique(stellar_parameters[name]) for name in names]))


    def neighbours(self, effective_temperature, surface_gravity, metallicity):
        """
        Return the indices of the neighbouring model points.
        """

        names = ("effective_temperature", "surface_gravity", "metallicity")
        values = (effective_temperature, surface_gravity, metallicity)
        nearest_upper_index = np.array([
            self._stellar_parameters[name].searchsorted(value) \
            for name, value in zip(names, values)])
        nearest_lower_index = np.array(nearest_upper_index) - 1

        nearest_upper_point = [self._stellar_parameters[name][i] \
            for name, i in zip(names, nearest_upper_index)]
        nearest_lower_point = [self._stellar_parameters[name][i] \
            for name, i in zip(names, nearest_lower_index)]

        indices = np.ones(self.stellar_parameters.size, dtype=bool)
        for name, lower, upper \
        in zip(names, nearest_lower_point, nearest_upper_point):
            indices *= (self.stellar_parameters[name] >= lower) \
                * (self.stellar_parameters[name] <= upper)

        num = indices.sum()
        if 8 > num:
            raise ValueError("not enough neighbouring points (8 > {}) to "
                "perform interpolation".format(num))

        elif num > 8:
            raise IndexError("too many neighbouring points ({} > 8) found"
                .format(num))

        return indices


    def __call__(self, *args, **kwargs):
        """ Alias to Interpolator.interpolate """
        return self.interpolate(*args, **kwargs)


    def interpolate(self, effective_temperature, surface_gravity, metallicity,
        opacity_scale="logTau5000"):
        """
        Return the interpolated photospheric quantities on a common opacity
        scale.
        """

        opacity_index = self.photosphere_cols.index(opacity_scale)

        try:
            indices = self.neighbours(effective_temperature, surface_gravity,
                metallicity)

        except (ValueError, IndexError):
            raise ValueError("cannot interpolate model atmosphere because grid"
                " is mal-formed")

        # Resample the opacities to a common scale on tau_5000 Angstroms
        photospheres = self.photospheres[indices].copy()
        opacities = common_opacity_scale(photospheres, opacity_index)

        # Put the stellar parameters on a unit cube scale
        stellar_parameters = self.stellar_parameters[indices]\
            .view(float).reshape(-1, 3)

        normed_subgrid = stellar_parameters - np.min(stellar_parameters, axis=0)
        p = np.array([effective_temperature, surface_gravity, metallicity]) \
            - np.min(stellar_parameters, axis=0)
        p /= np.max(normed_subgrid, axis=0)
        normed_subgrid /= np.max(normed_subgrid, axis=0)

        # Re-scale Pe and Pg to be logarithmic quantities?
        columns = [] + list(self.photosphere_cols)
        unlog_quantities = []
        logarithmic_quantities = ("Pe", "Pg")
        for quantity in logarithmic_quantities:
            try:
                index = columns.index(quantity)

            except ValueError:
                continue

            else:
                photospheres[:, :,  index] = np.log10(photospheres[:, :, index])
                columns[index] = "log{}".format(quantity)
                unlog_quantities.append((quantity, index))

        # Scale the radius?
        # [TODO]

        # Re-sample the photospheres onto a common opacity scale
        resampled_photospheres = np.zeros(photospheres.shape)
        for i, photosphere in enumerate(photospheres):
            resampled_photospheres[i] = \
                resample_photosphere(opacities, photosphere, opacity_index)

        # Then perform the interpolation using the empirical optimised
        # coefficients from Masseron (2006)
        interpolated_photosphere = np.zeros(self.photospheres.shape[1:])
        interpolated_photosphere[:, opacity_index] = opacities

        for j, quantity in enumerate(self.photosphere_cols):
            if j == opacity_index: continue

            # Scale using Masseron coefficients (where applicable)
            coefficients = _eval(self._coefficients.get(quantity, np.zeros(3)),
                {
                    "teff": effective_temperature,
                    "logg": surface_gravity,
                    "z": metallicity
                })
            scaled_range = np.abs(np.max(stellar_parameters, axis=0) \
                - np.min(stellar_parameters, axis=0)) / np.array([3200, 5, 4])
            scales = 1 - coefficients * scaled_range
            point = p.copy()**scales

            # Interpolate the quantities
            interpolated_photosphere[:, j] = scipy.interpolate.griddata(
                normed_subgrid, resampled_photospheres[:, :, j],
                point.reshape(1, 3)).flatten()

        # Rescale any logarithmic quantities
        for quantity, index in unlog_quantities:
            interpolated_photosphere[:, index] = \
                10**interpolated_photosphere[:, index]

        return np.core.records.fromarrays(interpolated_photosphere.T,
            names=self.photosphere_cols)




def resample_photosphere(opacities, photosphere, opacity_index):
    """ Resample photospheric quantities onto a new opacity scale. """

    resampled_photosphere = np.zeros(photosphere.shape)
    n_quantities = photosphere.shape[1]
    for i in range(n_quantities):
        if i == opacity_index: continue
        # Create spline function
        tck = scipy.interpolate.splrep(photosphere[:, opacity_index],
            photosphere[:, i])

        # Evaluate photospheric quantities at the new opacities
        resampled_photosphere[:, i] = scipy.interpolate.splev(opacities, tck)

    # Update photosphere with new opacities
    resampled_photosphere[:, opacity_index] = opacities
    return resampled_photosphere



def common_opacity_scale(photospheres, opacity_index):
    """ Returns a common opacity scale for the photospheres in question. """

    # Find the extent of the photospheres
    # model, depth_points, photospheric_properties
    pmin = photospheres[:, 0, opacity_index].max()
    pmax = photospheres[:, -1, opacity_index].min()

    # Pick the rescaled points by using the spacing information available in
    # the first model.
    first_model_scale = photospheres[0, :, opacity_index]
    fractional_shift = np.diff(first_model_scale)/np.ptp(first_model_scale)
    return np.hstack([pmin, pmin + np.cumsum(fractional_shift) * (pmax - pmin)])


def _eval_single(scale, env):
    if isinstance(scale, (int, float)):
        return scale
    default_env = { 
        "locals": None,
        "globals": None,
        "__name__": None,
        "__file__": None,
        "__builtins__": None,
    }
    default_env.update(env)
    return eval(scale, default_env)


def _eval(scales, env):
    return np.array([_eval_single(scale, env) for scale in scales])



if __name__ == "__main__":


    import matplotlib.pyplot as plt

    point = [5777., 4.445, 0.0]

    marcs = Interpolator()
    sun_interp = marcs.interpolate(*point)
    neighbours = marcs.neighbours(*point)

    fig, axes = plt.subplots(4)


    quantities = ("Depth", "T", "Pe", "Pg")
    for i, (ax, quantity) in enumerate(zip(axes, quantities)):
    
        quantity_ranges = marcs.photospheres[neighbours, :, i+1]
        min_range = np.min(quantity_ranges, axis=0)
        max_range = np.max(quantity_ranges, axis=0)

        ax.fill_between(sun_interp["logTau5000"], min_range, max_range,
            facecolor="#cccccc")
        ax.plot(sun_interp["logTau5000"], min_range, c="#666666")
        ax.plot(sun_interp["logTau5000"], max_range, c="#666666")
        
        ax.plot(sun_interp["logTau5000"], sun_interp[quantity], 'k', lw=2)

        ax.set_xlabel("tau(5000)")
        ax.set_ylabel(quantity)

        if i > 1:
            ax.set_yscale('log')

    print(sun_interp)

    plt.show()
    




