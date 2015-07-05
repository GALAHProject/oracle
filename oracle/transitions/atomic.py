#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" Represent atomic transitions. """

from __future__ import absolute_import, print_function

__all__ = ["AtomicTransition"]
__author__ = "Andy Casey <arc@ast.cam.ac.uk>"

import json
import yaml

from numpy import log, nan


class AtomicTransition(object):

    _DEFAULTS = {}

    def __init__(self, wavelength, species, e_low, log_gf, mask=None, **kwargs):
        """
        Initialise the class.
        """

        self.wavelength = float(wavelength)
        self.species = float(species)
        self.e_low = float(e_low)
        self.log_gf = float(log_gf)
        self.mask = mask

        self._raw = self._DEFAULTS.copy()
        self._raw.update({
            "wavelength": self.wavelength,
            "species": self.species,
            "e_low": self.e_low,
            "log_gf": self.log_gf,
            "mask": self.mask
        })
        self._raw.update(kwargs)

    def __str__(self):
        return unicode(self).encode("utf-8").strip()

    def __unicode__(self):
        return u"<{element} {ion} at {wavelength:.2f} Å>".format(
            element=self.element, ion="I"*self.ion, wavelength=self.wavelength)

    def __repr__(self):
        return "<oracle.transitions.{kls} {element} {ion} at {wavelength:.1f} "\
            "Angstroms at {location}>".format(kls=self.__class__.__name__,
                element=self.element, ion="I" * self.ion,
                wavelength=self.wavelength, location=hex(id(self)))

    def __eq__(self, transition):
        # If we have level information, that should be used.
        raise NotImplementedError

    def to_json(self, **kwargs):
        return json.dumps(self._raw, **kwargs)

    def to_yaml(self, **kwargs):
        return yaml.dump(self._raw, **kwargs)
    
    @property
    def element(self):
        periodictable = """H                                                  He
                           Li Be                               B  C  N  O  F  Ne
                           Na Mg                               Al Si P  S  Cl Ar
                           K  Ca Sc Ti V  Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr
                           Rb Sr Y  Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I  Xe
                           Cs Ba Lu Hf Ta W  Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn
                           Fr Ra Lr Rf Db Sg Bh Hs Mt Ds Rg Cn UUt"""
        
        lanthanoids    =  "La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb"
        actinoids      =  "Ac Th Pa U  Np Pu Am Cm Bk Cf Es Fm Md No"
        
        periodictable = periodictable.replace("Ba ", "Ba " + lanthanoids) \
            .replace("Ra ", "Ra " + actinoids).split()
        return periodictable[self.atomic_number - 1]

    @property
    def atomic_number(self):
        return int(self.species)

    @property
    def ion(self):
        return int((self.species % 1) * 10) + 1

    @property
    def equivalent_width(self):
        return self._raw.get("equivalent_width", nan)

    @equivalent_width.setter
    def equivalent_width(self, value):
        self._raw["equivalent_width"] = value

    @property
    def log_eps(self):
        return self._raw.get("log_eps", nan)

    @log_eps.setter
    def log_eps(self, value):
        self._raw["log_eps"] = value

    @property
    def reduced_equivalent_width(self):
        """
        Return the reduced equivalent width of the line, if the equivalent width
        has been measured.

        The reduced equivalent width is defined by:

        >> REW = log(equivalent width/wavelength)

        :raise AttributeError:
            If no equivalent width has been measured for this line.
        """
        return log(self.equivalent_width/self.wavelength)

    # Some default line behaviour that we don't want passed to _raw unless
    # explicitly defined.
    @property
    def continuum_degree(self):
        return self._raw.get("continuum_degree", 0) # Assumes normalised spectra

    @continuum_degree.setter
    def continuum_degree(self, value):
        self._raw["continuum_degree"] = value

    @property
    def fitting_region(self):
        return self._raw.get("fitting_region",
            [self.wavelength - 1.5, self.wavelength + 1.5])

    @fitting_region.setter
    def fitting_region(self, value):
        self._raw["fitting_region"] = value

    @property
    def v_rad_tolerance(self):
        return self._raw.get("v_rad_tolerance", 1.0) # km/s

    @v_rad_tolerance.setter
    def v_rad_tolerance(self, value):
        self._raw["v_rad_tolerance"] = value
