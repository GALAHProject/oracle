#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" Represent and model atomic transitions """

from __future__ import absolute_import, print_function

__all__ = ["AtomicTransition"]
__author__ = "Andy Casey <arc@ast.cam.ac.uk>"


import logging

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
            else "{0} {1} ".format(self.element, "I" * self.ionisation_level)
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
        constrain_parameters=None, atmosphere_kwargs=None,
        synthesis_kwargs=None, full_output=False, **kwargs):
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

        :param atmosphere_kwargs: [optional]
            Keyword arguments to supply to the atmosphere interpolator. This is
            only employed when a blending spectrum is required.

        :type atmosphere_kwargs:
            dict

        :param synthesis_kwargs: [optional]
            Keyword arguments to supply to the synthesis function. This is only
            employed when a blending spectrum is required.

        :type synthesis_kwargs:
            dict

        :param full_output:
            Return additional information about the resulting fit.

        :type full_output:
            bool

        :returns:
            The model parameters.

        :raises ValueError:
            - If an unavailable profile type was provided.
            - If `continuum_order` is not a positive integer.

        :raises TypeError:
            - If `continuum_order` is not an integer-like type.

        """

        if initial_theta is None:
            initial_theta = {}

        kind, available = kind.lower(), ("gaussian", "voigt", "lorentzian")
        if kind not in available:
            raise ValueError("available profile types are {}".format(available))

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

        if atmosphere_kwargs is None:
            atmosphere_kwargs = {}

        if synthesis_kwargs is None:
            synthesis_kwargs = {}

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
                for i in range(continuum_order + 1)]

        # Blending spectrum FWHM?
        if self.blending_transitions is not None \
        and len(self.blending_transitions) > 0:
            parameters.append("blending_fwhm")

            # Since we have to synthesise a blending spectrum, check that we
            # have the required stellar parameters in initial_theta
            stellar_parameters = ("effective_temperature", "surface_gravity",
                "metallicity", "microturbulence")
            if not all([k in initial_theta for k in stellar_parameters]):
                raise KeyError("stellar parameters ({}) are required for "\
                    "synthesis because this transition has blending transitions"
                    .format(", ".join(stellar_parameter_keys)))

        logger.debug("{0} has {1} parameters: {2}".format(self, len(parameters),
            ", ".join(parameters)))
            

        # Synthesise a blending spectrum if necessary

        # Estimate the initial values of parameters that were not given.

        # Create a copy of the data and (if it exists) apply a mask to the data.


        # Define the minimisation function
            # -> check for any violations in constrain_parameters or sensible
            #    physics
            # -> copy the blending spectrum
            # -> smooth the blending spectrum
            # -> form the absorption profile (from unity)
            # -> multiply the two together
            # -> scale it by the continuum

            # -> [eventually] account for outlier pixels
            # -> calculate a chi_sq value


        """
        Simplest case:

        transition = AtomicTransition(wavelength=5775.12)
        transition.fit(data)

        # What do
        # No continuum
        # No synthesis
        # No mask
        # No species information
        # ... Fit with a Gaussian profile.
        #     Return the fit parameters.
        #     Integrate the profile for an equivalent width, which we store in
        #        the class.


        # With species information:

        transition = AtomicTransition(wavelength=5775.12, excitation_potential=1.2,
            species=26.0, loggf=0.00)
        transition.fit(data)

        # What do?
        # Same as above.
        # We can't provide an abundance unless we have stellar parameter
        # information as well.

        # Perhaps self.abundance should be self.abundance(effective_temperature,
            surface_gravity, metallicity, ...etc)?
        # --> Yes.
    

        # With a mask:

        transition = AtomicTransition(wavelength=5775.12, excitation_potential=1.2,
            species=26.0, loggf=0.00, mask=[5776.5, 5779.1])
        transition.fit(data)

        # What do?
        # Same as above, except use the mask in the fitting process


        # With blending transitions:

        transition = AtomicTransition(wavelength=5775.12, excitation_potential=1.2,
            species=26.0, loggf=0.00, mask=[5776.5, 5779.1],
            blending_transitions=[[5774.1, 12.0, 0.00, 0.00]])
        
        [or AtomicTransitions could be given as the blending transitions]
        transition.fit(data)

        # raise an error: stellar parameters are required


        # With blending transitions but no atomic data for the line itself?
        # --> That's OK, but .abundance() will fail


        # OK, so we have blending transitions. Now let's see what we do for
        # fit_profile()
        # options:
        # kind = gaussian
        # continuum_order = None,
        # constrain_parameters={'wavelength'}
        # atmosphere kwargs
        # outlier_treatment

        # free parameters are:
        # line depth, fwhm, fwhm of blending spectrum, any continuum parameters,
        """




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
        ionisation_level = int(10 * (species - atomic_number))
        
    return (element, atomic_number, ionisation_level, species)
        
