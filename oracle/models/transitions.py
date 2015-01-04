#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" Represent and model atomic transitions """

from __future__ import absolute_import, print_function

__all__ = ["AtomicTransition"]
__author__ = "Andy Casey <arc@ast.cam.ac.uk>"

import logging

import numpy as np
from scipy import ndimage, optimize as op

from oracle import synthesis
from oracle.models import profiles
from oracle.specutils import Spectrum1D
from oracle.utils import (atomic_number as parse_atomic_number,
    element as parse_element)

logger = logging.getLogger("oracle")


class AtomicTransition(object):

    def __init__(self, wavelength, species=None, excitation_potential=None,
        loggf=None, blending_transitions=None, mask=None, **kwargs):
        """
        Initialise the class
        """

        self.wavelength = float(wavelength)

        float_if_not_none = lambda x: float(x) if x is not None else x
        self.blending_transitions = blending_transitions
        if blending_transitions is not None:
            self.blending_transitions = np.atleast_2d(blending_transitions)
        self.excitation_potential = float_if_not_none(excitation_potential)
        self.loggf = float_if_not_none(loggf)
        self.mask = mask

        # Species can be 'Fe I', '26.0', 26.1, 'Fe' (assumed neutral)
        self.element, self.atomic_number, self.ionisation_level, self.species \
            = _parse_species(species)

        self._init_kwargs = kwargs.copy()


    def __str__(self):
        return unicode(self).encode("utf-8").strip()


    def __unicode__(self):
        species_repr = "Unspecified transition" if self.species is None \
            else "{0} {1}".format(self.element, "I" * self.ionisation_level)
        return u"{species_repr} at {wavelength:.1f} Å".format(
            species_repr=species_repr, wavelength=self.wavelength)


    def __repr__(self):
        species_repr = "" if self.species is None else "{0} {1} ".format(
            self.element, "I" * self.ionisation_level)
        # [TODO] astropy.units to show wavelength units
        return "<{module}.AtomicTransition {species_repr}at {wavelength:.1f} "\
            "Angstrom at {location}>".format(module=self.__module__,
                species_repr=species_repr, wavelength=self.wavelength,
                location=hex(id(self)))


    def fit_profile(self, data, initial_theta=None, kind="gaussian", 
        continuum_order=None, surrounding=1.5, outlier_pixels=None,
        constrain_parameters=None, synthesise_kwargs=None, optimise_kwargs=None,
        full_output=False, **kwargs):
        """
        Model and fit the atomic transition with an absorption profile. This
        function can account for continuum, radial velocity offsets (by limiting
        the wavelength range), blending transitions, and outlier pixels.

        While this approach can account for blended lines by synthesising them,
        the user must remember that opacities are not strictly employed properly
        using this method.

        :param data:
            The observed data. This spectrum should include the wavelength of
            the transition.

        :type data:
            :class:`oracle.specutils.Spectrum1D`

        :param initial_theta: [optional]
            Initial estimates of the model parameters :math:`\Theta`. The number
            of model parameters varies depending on the optional inputs, but the
            only common parameter is `line_depth`. The wavelength can
            be a parameter (`wavelength`) if `constrain_parameters` has
            some constraint for `wavelength`, otherwise the wavelength is fixed
            and is not a model parameter.

            The profile parameters will vary depending on the `kind`:

                Gaussian: `fwhm`
                Voigt: `fwhm`, `shape`
                Lorentzian: `scale`

            Where applicable, continuum parameters are represented by
            `continuum.0`, ..., `continuum.N`. If a blending spectrum is to be
            produced (e.g., when `blending_transitions` is not None) then the
            FWHM of the Gaussian kernel to smooth the blending spectrum is
            also a parameter, accessible through `blending_fwhm`.

            :note:
            If `blending_transitions` were supplied to the `AtomicTransition`
            class, then a spectrum will be synthesised. In this situation the
            stellar parameters of the parent model must be provided. These keys
            are `effective_temperature`, `surface_gravity`, `metallicity`, and
            `microturbulence`.

            # Discuss outlier parameters.

        :type initial_theta:
            dict

        :param kind: [optional]
            The kind of profile function to use. Available options are Gaussian,
            Voigt or Lorentzian.

        :type kind:
            str

        :param continuum_order: [optional]
            The order of the polynomial to use to fit the continuum. By default,
            this option is set to None, meaning that no continuum is included.
            In this case the data are assumed to be continuum-normalised (e.g.,
            the absorption begins at unity). Setting this value to zero implies
            one continuum parameter a scaling offset, setting to one implies two
            parameters (ax + b),and so on.

            If continuum is treated, initial values for the continuum can be
            provided by the 'continuum.0', ..., 'continuum.N' keys in
            `initial_theta`.

        :type continuum_order:
            int

        :param surrounding: [optional]
            The region of spectrum to model either side of the `wavelength`. By
            default this is 1.5 Angstroms.

        :type surrounding:
            float

        :param outlier_pixels: [optional]
            Decide how to model the outlier pixels. Currently: nothing.

        :type outlier_pixels:
            str

        :param constrain_parameters: [optional]
            Specify constraints for the model parameters. Model parameters are
            expected as keys, and their constraints (lower, upper) are expected
            as values. Use `None` in the lower or upper bound to specify no
            constraint on that edge.

        :type constrain_parameters:
            dict

        :param synthesise_kwargs: [optional]
            Keyword arguments to supply to the synthesise function. This is only
            employed when a blending spectrum is required.

        :type synthesise_kwargs:
            dict

        :param optimise_kwargs: [optional]
            Keyword arguments to supply to the optimisation function.

        :type optimise_kwargs:
            dict

        :param full_output:
            Return additional information about the resulting fit.

        :type full_output:
            bool

        :returns:
            The model parameters.
        """

        if not isinstance(data, Spectrum1D):
            raise TypeError("data must be a oracle.specutils.Spectrum1D object")

        if not (data.disp[-1] > self.wavelength > data.disp[0]):
            raise ValueError("the atomic transition is not in the bounds of the"
                " input data ({0:.0f} outside [{1:.0f}, {2:.0f}])".format(
                    self.wavelength, data.disp[0], data.disp[-1]))

        if initial_theta is None:
            initial_theta = {}

        kind, available = kind.lower(), ("gaussian", "voigt", "lorentzian")
        if kind not in available:
            raise ValueError("available profile types are {}".format(available))

        # [TODO] fix this
        if kind in ("voigt", "lorentzian"):
            raise NotImplementedError("sry gooby")

        if continuum_order is not None:
            try:
                continuum_order = int(continuum_order)
            except (ValueError, TypeError):
                raise TypeError("continuum order must be an integer")
            if continuum_order < 0:
                raise ValueError("continuum order must be a positive integer")

        try:
            surrounding = float(surrounding)
        except (ValueError, TypeError):
            raise TypeError("surrounding region must be a float")
        if surrounding < 0:
            raise ValueError("surrounding region must be a positive float")

        # outlier_pixels

        if constrain_parameters is None:
            constrain_parameters = {}
        else:
            # Set lower case keys.
            constrain_parameters = dict((k.lower(), v) \
                for k, v in constrain_parameters.iteritems())

        if synthesise_kwargs is None:
            synthesise_kwargs = {}

        # Establish the fit parameters
        parameters = ["line_depth"]
        profile_parameters = {
            "gaussian": ["fwhm"],
            "lorentzian": ["scale"],
            "voigt": ["fwhm", "shape"]
        }
        parameters.extend(profile_parameters[kind])

        # Is wavelength a free-ish parameter, or completely fixed?
        if "wavelength" in constrain_parameters:
            parameters.append("wavelength")

        # Continuum
        if continuum_order is not None:
            parameters.extend(["continuum.{}".format(i) \
                for i in range(continuum_order + 1)])

        # Blending spectrum FWHM?
        stellar_parameter_keys = ("effective_temperature", "surface_gravity",
            "metallicity", "microturbulence")
        if self.blending_transitions is not None \
        and len(self.blending_transitions) > 0:
            parameters.append("blending_fwhm")

            # Since we have to synthesise a blending spectrum, check that we
            # have the required stellar parameters in initial_theta
            if not all([k in initial_theta for k in stellar_parameter_keys]):
                raise KeyError("stellar parameters ({}) are required for "\
                    "synthesis because this transition has blending transitions"
                    .format(", ".join(stellar_parameter_keys)))

        logger.debug("{0} has {1} parameters: {2}".format(self, len(parameters),
            ", ".join(parameters)))

        # Create a copy of a region of the data and (where applicable,..) apply
        # a mask to the data.
        region = (self.wavelength - surrounding, self.wavelength + surrounding)
        data = data.slice(region, copy=True)
        if self.mask is not None:
            data = data.mask_dispersions(self.mask)
            
        # Estimate the initial values of parameters that were not given.
        missing_parameters = set(parameters).difference(initial_theta)
        logger.debug("Need to estimate {0} parameters: {1}".format(
            len(missing_parameters), ", ".join(missing_parameters)))

        # Order of estimates (in increasing difficulty)
        # -> wavelength, continuum.N, line_depth, fwhm (et al.), blending_fwhm
        theta = { "wavelength": self.wavelength }
        for parameter in set(parameters).intersection(initial_theta):
            # Use setdefault so we don't overwrite the wavelength
            theta.setdefault(parameter, initial_theta[parameter])

        # Synthesise a blending spectrum if necessary
        if "blending_fwhm" in parameters:
            effective_temperature, surface_gravity, metallicity, microturbulence\
                = [initial_theta[k] for k in stellar_parameter_keys]

            pixel_size = np.diff(data.disp).min()
            synthesise_kwargs.setdefault("oversample", 4)
            synthesise_kwargs.setdefault("wavelength_step", pixel_size)

            blending_spectrum = Spectrum1D(*synthesis.synthesise(
                effective_temperature, surface_gravity, metallicity,
                microturbulence, self.blending_transitions, region[0], region[1],
                **synthesise_kwargs))
        else:
            blending_spectrum = None

        # Do we need to estimate continuum parameters?
        # Note here it's either all or nothing. Either initial_theta provides us
        # with all of the continuum.* parameters we need, or we estimate them
        # all ourselves
        if any([p.startswith("continuum.") for p in missing_parameters]):
            if blending_spectrum is None:
                # Polyfit
                # [TODO] any outlier removal required?
                coefficients = np.polyfit(data.disp, data.flux, continuum_order)

            else:
                # Divide with the blending spectrum
                disp_matrix = np.ones((continuum_order+1, data.disp.size))
                for i in range(continuum_order + 1):
                    disp_matrix[i] *= data.disp**i

                # Ensure blending spectrum is on the same pixel scale (e.g.,
                # there has been no rebinning)
                rbs = np.interp(
                    data.disp, blending_spectrum.flux, blending_spectrum.flux)
                A = (data.flux/rbs).T
                coefficients = np.linalg.lstsq(disp_matrix.T, A)[0][::-1]

            # Update theta with the coefficients we have estimated
            theta.update(dict(zip(["continuum.{}".format(i) \
                for i in range(continuum_order + 1)], coefficients)))

        else:
            # Get the continuum coefficients, we may need them later.
            if continuum_order is None:
                coefficients = [1] # sneaky sneaky

            else:
                coefficients = [theta["continuum.{}".format(i)] \
                    for i in range(continuum_order + 1)]

        # Do we need to estimate the line depth?
        if "line_depth" in missing_parameters:
            continuum_at_wavelength = np.polyval(coefficients, self.wavelength)
            index = data.disp.searchsorted(self.wavelength)
            line_depth = 1 - data.flux[index]/continuum_at_wavelength
            
            # Ensure the initial value is limited within a sensible range
            tolerance = kwargs.pop("initial_line_depth_tolerance", 1e-3)
            theta["line_depth"] = np.clip(line_depth, tolerance, 1 - tolerance)

        # Do we need to estimate any of the profile parameters?
        # [TODO] here we just do Gaussian because simplicity. Expand on that
        if "fwhm" in missing_parameters:
            mid_flux = continuum_at_wavelength * (1 - 0.5 * theta["line_depth"])

            # Find the wavelengths from the central wavelength where this value
            # exists
            pos_mid_point, neg_mid_point = np.nan, np.nan
            index = data.disp.searchsorted(self.wavelength)
            pos_indices = data.flux[index:] > mid_flux
            if pos_indices.any():
                pos_mid_point = data.disp[index:][pos_indices][0]

            neg_indices = data.flux[:index] > mid_flux
            if neg_indices.any():
                neg_mid_point = data.disp[:index][neg_indices][-1]

            # If both positive and negative mid points are nans, then we cannot
            # estimate the FWHM with this method
            if not np.isfinite([pos_mid_point, neg_mid_point]).any():
                raise ValueError("cannot estimate initial fwhm")

            # If either mid point is a nan, take the initial fwhm from one side
            # of the profile
            initial_fwhm = pos_mid_point - neg_mid_point
            if not np.isfinite(initial_fwhm):
                initial_fwhm = np.array([
                    2 * (pos_mid_point - self.wavelength),
                    2 * (self.wavelength - neg_mid_point)
                ])
                finite = np.isfinite(initial_fwhm)
                initial_fwhm = initial_fwhm[finite][0]

            theta["fwhm"] = initial_fwhm

        # Do we need to estimate the blending FWHM? If so we will just take it
        # from the initial profile FWHM
        if "blending_fwhm" in missing_parameters:
            theta["blending_fwhm"] = theta["fwhm"]

        assert len(set(parameters).difference(theta)) == 0, "Missing parameter!"

        invalid_return = np.nan

        def chi_sq(x, full_output=False):

            # Extra variables for free due to scope:
            # invalid_return, constrain_parameters, blending_spectrum, data,
            # continuum_order, Spectrum1D?

            xd = dict(zip(parameters, x))

            # Physically sensible constraints
            if 0 > xd["fwhm"] or not (1 > xd["line_depth"] > 0) \
            or 0 > xd.get("blending_fwhm", 1):
                return invalid_return \
                    if not full_output else 3 * (invalid_return, )

            # Constraints provided by the user
            for parameter, (lower, upper) in constrain_parameters.iteritems():
                if (lower is not None and lower > xd[parameter]) \
                or (upper is not None and upper < xd[parameter]):
                    return invalid_return \
                        if not full_output else 3 * (invalid_return, )

            # Smooth the blending spectrum
            if blending_spectrum is not None:
                fwhm = xd["blending_fwhm"]/np.diff(blending_spectrum.disp).min()
                blending_spectrum_flux = ndimage.gaussian_filter(
                    blending_spectrum.flux, fwhm/2.355)

                # Resample to the data dispersion points
                blending_spectrum_flux = np.interp(data.disp,
                    blending_spectrum.disp, blending_spectrum_flux)
            else:
                blending_spectrum_flux = np.ones(data.disp.size)

            # Form the absorption profile
            absorption_profile = \
                1 - xd["line_depth"] * profiles.gaussian(
                    xd.get("wavelength", self.wavelength), xd["fwhm"]/2.355,
                    data.disp)

            # Scale by the continuum (if applicable)
            continuum = 1 if continuum_order is None \
                else np.polyval([xd["continuum.{}".format(i)] \
                    for i in range(continuum_order + 1)], data.disp)

            # Put everything together 
            model = blending_spectrum_flux * absorption_profile * continuum

            # Calculate the chi-squared value
            chi_sqs = (model - data.flux)**2 * data.ivariance
            finite = np.isfinite(chi_sqs)
            chi_sq = chi_sqs[finite].sum()

            if full_output:
                r_chi_sq = chi_sq / (finite.sum() - len(xd) - 1)
                return (chi_sq, r_chi_sq, Spectrum1D(data.disp, model))
            return chi_sq

        minimisation_function = chi_sq

        # Check for any violations in the initial parameters
        p0 = np.array([theta[p] for p in parameters])
        try:
            # initial chi-sq, initial reduced chi-sq, initial model fit
            ic, irc, imf = minimisation_function(p0, full_output=True)

        except:
            logger.exception("Exception raised when trying to start profile "
                "fitting:")
            raise

        else:
            if not np.isfinite(ic) or ic == invalid_return:
                raise ValueError("invalid value returned by minimisation "
                    "function when initial theta values were provided")

        # [TODO] Allow for something other than ye-olde Nelder & Meade?
        op_info = { "p0": theta }
        p1, fopt, num_iter, num_funcalls, warnflag = op.fmin(
            minimisation_function, p0, disp=False, full_output=True)
        oc, orc, omf = minimisation_function(p1, full_output=True)
        op_info.update({
            "fopt": fopt,
            "num_iter": num_iter,
            "num_funcalls": num_funcalls,
            "warnflag": warnflag
        })
        optimal_theta = dict(zip(parameters, p1))

        # Any sanity checks?
        if "blending_fwhm" in optimal_theta:
            fwhm_change = optimal_theta["blending_fwhm"]/optimal_theta["fwhm"]
            if not (2 > fwhm_change > 0.5):
                logging.debug("Profile and blending FWHM differ substantially "
                    "from each other: profile_fwhm = {0:.2f} Angstrom, blending"
                    "_fwhm = {1:.2f}".format(optimal_theta["fwhm"],
                        optimal_theta["blending_fwhm"]))

        # What information should we save to the class?
        self._fit_profile_result = (optimal_theta, oc, orc, omf, op_info)

        profile_integrals = {
            # Other profile functions will probably need additional information
            # other than just 'x' (the optimal dictionary)
            # For those playing at home: pi^2/2.355 == 0.752634331594699
            "gaussian": lambda x: x["line_depth"] * abs(x["fwhm"]) * 0.752634332
        }
        # Calculate the equivalent width (in milli Angstroms)
        # [TODO] use astropy.units to provide some relative measure
        self.equivalent_width = 1000 * profile_integrals[kind](optimal_theta)

        if full_output:
            return (optimal_theta, oc, orc, omf, op_info)

        # I could not sleep for 65 days because I worried about whether I should
        # be returning a list of values or a dictionary by default. In the end I
        # decided that since the length (and order) of model parameters can vary
        # depending on the input keywords, and that the model parameters are not
        # easily accessible, I should return a dictionary here to be 'explicit'.

        # (The above paragraph pleases my OCD)

        return optimal_theta


def _parse_ionisation_level(ionisation_level):
    """
    Common imputs:
    1, 2, I, II
    """

    if isinstance(ionisation_level, int):
        return ionisation_level

    # String assumed
    ionisation_level = ionisation_level.upper()
    if "I" in ionisation_level:
        # Assuming less than ionisation state 4 (e.g., IV)
        return ionisation_level.count("I")
    return int(ionisation_level)


def _parse_species(species):

    if species is None:
        # No idea.
        return (None, None, None, None)

    try:
        # Float species are easiest
        species = float(species)

    except (ValueError, TypeError):
        # Element string representations are the only thing that we will accept
        # now

        atomic_number = parse_atomic_number(species)
        element = parse_element(atomic_number)

        # Unless provided, ground state is assumed.
        ionisation_level = _parse_ionisation_level(species.split()[1]) \
            if " " in species else 1

        species = atomic_number + ionisation_level/10.
        return (element, atomic_number, ionisation_level, species)


    else:
        element = parse_element(species)
        atomic_number = int(species)
        ionisation_level = int(10 * (species - atomic_number)) + 1
        
    return (element, atomic_number, ionisation_level, species)
        
