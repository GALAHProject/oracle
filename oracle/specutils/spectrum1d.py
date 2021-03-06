#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" Class for dealing with 1D spectra """

from __future__ import division, print_function

__all__ = ["Spectrum1D", "Spectrum"]
__author__ = "Andy Casey <arc@ast.cam.ac.uk>"

# Standard library
import logging
import os
import warnings
from pkg_resources import resource_stream

# Third-party
import numpy as np
import yaml
import astropy.units as u
from astropy.io import fits
from astropy.time import Time

# This module.
from . import helcorr

logger = logging.getLogger("oracle")

# Warn about missing variance arrays, but only once.
class MissingVarianceWarning(Warning):
    pass
warnings.simplefilter("once", MissingVarianceWarning)

class Spectrum1D(object):
    """
    This is a temporary class holder for a Spectrum1D object until the
    :class:`astropy.specutils.Spectrum1D` module has advanced sufficiently to
    replace it.
    """
    
    def __init__(self, disp, flux, variance=None, headers=None, **kwargs):
        """
        Initializes a `Spectrum1D` object with the given dispersion, flux, and
        variance arrays.

        :param disp:
            The dispersion points (e.g., wavelengths) of the spectrum.

        :type disp:
            :class:`np.array`

        :param flux:
            The associated flux points of the spectrum.

        :type flux:
            :class:`np.array`

        :param variance: [optional]
            The variance of the flux points. If not provided, the variance in
            the flux is assumed to be Poisson-distributed.

        :type variance:
            :class:`np.array`

        :param headers: [optional]
            The metadata associated with the spectrum.

        :type headers:
            dict
        """

        if len(disp) != len(flux):
            raise ValueError("dispersion and flux must have the same length")

        if len(disp) == 0:
            raise ValueError("dispersion and flux cannot be empty arrays")
        
        self.disp = disp.copy()
        self.flux = flux.copy()
        if variance is None:
            warnings.warn("No variance array provided. Unless otherwise given,"\
                " the noise in spectra are assumed to be Poisson-distributed.",
                MissingVarianceWarning)
            self.variance = self.flux.copy()

        else:
            self.variance = variance.copy()

        # Zero or negative fluxes should have a high variance
        self.variance[0 >= self.flux] = kwargs.pop("zero_flux_variance", 1e6)

        # Minimum variance
        self.variance[0 >= self.variance] = kwargs.pop("min_variance", 1e-6)

        self.ivariance = 1.0/self.variance
        if headers is not None:
            self.headers = headers
        else:
            self.headers = {}
        return None


    def copy(self):
        """ Creates a copy of the object """
        return self.__class__(self.disp.copy(), self.flux.copy(),
            variance=self.variance.copy(), headers=self.headers.copy())
    

    def __contains__(self, value):
        try:
            # for AtomicTransitions
            return self.disp[-1] >= value.wavelength >= self.disp[0]
        except AttributeError:
            return self.disp[-1] >= value >= self.disp[0]


    @property
    def v_helio(self):
        """
        Return the heliocentric velocity, given the information available in
        the headers.
        """

        if hasattr(self, "_v_helio"):
            return self._v_helio

        # Needs: RA, DEC (degrees, J2000)
        # MJD of middle of exposure
        # Some observatory information

        # ORIGIN: Observatory
        #RA
        # DEC
        #UTSTART = '14:19:02'           / UT start                                       
        #UTEND   = '14:39:03'           / UT end      

        #ALT_OBS =                 1164 / Altitude of observatory in metres              
        #LAT_OBS =            -31.27704 / Observatory latitude in degrees                
        #LONG_OBS=             149.0661 / Observatory longitude in degrees  

        # MEANRA (degrees)
        # MEANDEC (degrees)

        # Get the observatory.
        alt_obs = self.headers.get("ALT_OBS", None)
        lat_obs = self.headers.get("LAT_OBS", None)
        long_obs = self.headers.get("LONG_OBS", None)
        if None in (alt_obs, lat_obs, long_obs):
            # Try and determine it from the observatory name, if it exists.
            origin = self.headers.get("ORIGIN", None)

            if origin is None:
                raise KeyError("no observatory information available (ALT_OBS, "
                    "LAT_OBS, LONG_OBS) or ORIGIN")

            with resource_stream(__name__, "observatories.yaml") as fp:
                observatories_dictionary = yaml.load(fp)

            origin = origin.strip().lower()
            if origin not in observatories_dictionary:
                raise KeyError("could not find {} in the observatory dictionary"\
                    .format(origin))

            observatory = observatories_dictionary[origin]
            alt_obs = observatory["altitude"]
            lat_obs = observatory["latitude"]
            raise WTFError()

        # Get the RA/DEC.
        ra = self.headers.get("RA", None) # assuming degrees
        dec = self.headers.get("DEC", None)
        if None in (ra, dec):
            ra = self.headers.get("MEANRA", None)
            dec = self.headers.get("MEANDEC", None)

        if None in (ra, dec):
            raise KeyError("no position information (looked for RA/DEC, MEANRA"\
                "/MEANDEC")

        # Time of observation.
        for k in ("UTDATE", "UTSTART", "UTEND"):
            if k not in self.headers:
                raise KeyError("cannot find key {} in headers".format(k))

        ut_start = Time("{0}T{1}".format(self.headers["UTDATE"].replace(":", "-"),
            self.headers["UTSTART"]), format="isot", scale="utc")
        ut_end = Time("{0}T{1}".format(self.headers["UTDATE"].replace(":", "-"),
            self.headers["UTEND"]), format="isot", scale="utc")

        # Get the MJD of the mid-point of the observation.
        mjd = (ut_end - ut_start).jd/2 + ut_start.mjd

        # Calculate the correction.
        v_helio = helcorr.helcorr(long_obs, lat_obs, alt_obs, ra, dec, mjd)
        self._v_helio = v_helio
        return v_helio



    def mask_by_dispersion(self, mask, mask_value=np.nan):

        flux = self.flux.copy()
        for start, end in mask:
            indices = self.disp.searchsorted([start, end])
            flux[indices[0]:indices[1]] = mask_value

        return self.__class__(self.disp, flux, variance=self.variance, headers=self.headers.copy())

    def slice(self, wavelengths, copy=False):
        """
        Slice a spectrum by some wavelengths.
        """

        assert len(wavelengths) == 2

        wavelength_start, wavelength_end = wavelengths
        index_start = self.disp.searchsorted(wavelength_start, side="left")
        index_end = self.disp.searchsorted(wavelength_end, side="right")

        if copy:
            disp = self.disp[index_start:index_end].copy()
            flux = self.flux[index_start:index_end].copy()
            variance = self.variance[index_start:index_end].copy()

        else:
            disp = self.disp[index_start:index_end]
            flux = self.flux[index_start:index_end]
            variance = self.variance[index_start:index_end]

        return self.__class__(disp, flux, variance=variance,
            headers=self.headers.copy())


    @classmethod
    def load_GALAH(cls, filename, normalised=False, rest=True, clean_edges=True,
        **kwargs):
        """
        Load a Spectrum1D object from the GALAH standardised FITS format.

        :param filename:
            The path of the filename to load.

        :type filename:
            str
        """

        clean_edges_limit = kwargs.pop("clean_edges_limit", 10)
        
        image = fits.open(filename, **kwargs)

        data_ext = 0 if not normalised else 2

        flux = image[data_ext].data
        if flux.size == 0 and normalised:
            raise ValueError("no normalised spectrum found")

        variance = image[data_ext + 1].data

        disp = image[data_ext].header["CRVAL1"] \
            + (np.arange(flux.size) - image[data_ext].header.get("CRPIX1", 0)) \
            * image[data_ext].header["CDELT1"]

        header_columns = ["CCD", "WG6_HASH", "NAME", "RA", "DEC",
            "PMRA", "PMDEC", "MAG", "DESCR", "FIBRE", "MOON_DEG", "V_HELIO"]

        headers = dict(zip(header_columns, [image[0].header.get(k, None) \
            for k in header_columns]))

        if clean_edges:
            # Look for really sharp changes at the edges of the spectrum
            y = np.abs(np.diff(flux))
            stds = (y - np.median(y))/np.std(y)

            lhs = np.any(stds[:50] > clean_edges_limit)
            if lhs:
                lhs_index = np.argmax(stds[:50]) + 1
            else:
                lhs_index = 0

            rhs = np.any(stds[-50:] > clean_edges_limit)
            if rhs:
                rhs_index = flux.size - np.argmax(stds[-50::][::-1]) - 1
            else:
                rhs_index = None

            disp = disp[lhs_index:rhs_index]
            flux = flux[lhs_index:rhs_index]
            if variance is not None:
                variance = variance[lhs_index:rhs_index]

        if rest:
            vrad = image[4].header.get("VRAD", np.nan)
            if np.isfinite(vrad):
                disp *= (1 - vrad/299792.458)
            else:
                logger.warn("No velocity information found; spectrum may not be at rest!")

        return cls(disp, flux, variance=variance, headers=headers)


    @classmethod
    def load(cls, filename, **kwargs):
        """Load a Spectrum1D from a given filename.
        
        :param filename:
            The path of the filename to load. Can be either simple 1D FITS files
            or an ASCII filename.

        :type filename:
            str
        """
        
        if not os.path.exists(filename):
            raise IOError("filename {0} does not exist" .format(filename))
        
        if filename.endswith('.fits'):
            image = fits.open(filename, **kwargs)
            
            header = image[0].header
            
            # Check for a tabular data structure
            if len(image) > 1 and image[0].data is None:

                names = [name.lower() for name in image[1].data.names]
                dispersion_key = 'wave' if 'wave' in names else 'disp'
                
                disp, flux = image[1].data[dispersion_key], image[1].data['flux']

                if 'error' in names or 'variance' in names:
                    variance_key = 'error' if 'error' in names else 'variance'
                    variance = image[1].data[variance_key]

            else:

                # According to http://iraf.net/irafdocs/specwcs.php ....
                #li = a.headers['LTM1_1'] * np.arange(a.headers['NAXIS1']) \
                #       + a.headers['LTV1']
                #a.headers['CRVAL1'] + a.headers['CD1_1'] * (li - 
                #a.headers['CRPIX1'])

                if np.all([key in header.keys() \
                    for key in ('CDELT1', 'NAXIS1', 'CRVAL1')]):

                    li = header.get("LTM1_1", 1) * np.arange(header['NAXIS1']) \
                        - header.get("LTV1", 0)

                    disp = header['CRVAL1']\
                        + header['CDELT1'] * (li - header.get("CRPIX1", 0))
            
                #if "LTV1" in header.keys():
                #    disp -= header['LTV1'] * header['CDELT1']

                #disp -= header['LTV1'] if header.has_key('LTV1') else 0
                flux = image[0].data
                variance = None
            
                # Check for logarithmic dispersion
                if "CTYPE1" in header.keys() and header["CTYPE1"] == "AWAV-LOG":
                    disp = np.exp(disp)

            # Add the headers in
            headers = {}
            for row in header.items():
                key, value = row
                
                # Check the value is valid
                try:
                    str(value)

                except TypeError:
                    logger.debug("Skipping header key {0}".format(key))
                    continue

                if len(key) == 0 or len(str(value)) == 0: continue
    
                if key in headers.keys():
                    if not isinstance(headers[key], list):
                        headers[key] = [headers[key]]
                    
                    headers[key].append(value)

                else:
                    headers[key] = value

            for key, value in headers.iteritems():
                if isinstance(value, list):
                    headers[key] = "\n".join(map(str, value))

        else:
            headers = {}
            loadtxt_kwargs = kwargs.copy()
            loadtxt_kwargs["unpack"] = True
            try:
                disp, flux, variance = np.loadtxt(filename, **loadtxt_kwargs)
            except:
                disp, flux = np.loadtxt(filename, **loadtxt_kwargs)
            
        # Clean the edges?
        clean_edges = kwargs.pop("clean_edges", True)
        if clean_edges:
            # Look for really sharp changes at the edges of the spectrum
            y = np.abs(np.diff(flux))
            stds = (y - np.median(y))/np.std(y)

            lhs = np.any(stds[:50] > 10)
            if lhs:
                lhs_index = np.argmax(stds[:50]) + 1
            else:
                lhs_index = 0

            rhs = np.any(stds[-50:] > 10)
            if rhs:
                rhs_index = flux.size - np.argmax(stds[-50::][::-1]) - 1
            else:
                rhs_index = None

            disp = disp[lhs_index:rhs_index]
            flux = flux[lhs_index:rhs_index]
            if variance is not None:
                variance = variance[lhs_index:rhs_index]

        return cls(disp, flux, variance=variance, headers=headers)


    def save(self, filename, clobber=True, **kwargs):
        """
        Save the `Spectrum1D` object to the specified filename.
        
        :param filename:
            The filename to save the Spectrum1D object to.

        :type filename:
            str

        :param clobber: [optional]
            Whether to overwrite the `filename` if it already exists.

        :type clobber:
            bool

        :raises IOError:
            If the filename exists and we were not asked to clobber it.
        """
        
        if os.path.exists(filename) and not clobber:
            raise IOError("filename '{0}' exists and we have been asked not to"\
                " clobber it".format(filename))
        
        if not filename.endswith('fits'):
            # ASCII
            data = np.hstack([
                self.disp.reshape(-1, 1),
                self.flux.reshape(-1, 1),
                self.variance.reshape(-1, 1)
                ])
            return np.savetxt(filename, data, **kwargs)
            
        else:          
            # Create a tabular FITS format
            disp = fits.Column(name='disp', format='1D', array=self.disp)
            flux = fits.Column(name='flux', format='1D', array=self.flux)
            var = fits.Column(name='variance', format='1D',
                array=self.variance)
            table_hdu = fits.new_table([disp, flux, var])

            # Create Primary HDU
            hdu = fits.PrimaryHDU()

            # Update primary HDU with headers
            for key, value in self.headers.iteritems():
                if len(key) > 8: # To deal with ESO compatibility
                    hdu.header.update('HIERARCH {}'.format(key), value)

                try:
                    hdu.header.update(key, value)
                except ValueError:
                    logger.warn("Could not save header key/value combination: "\
                        "{0} = {1}".format(key, value))

            # Create HDU list with our tables
            hdulist = fits.HDUList([hdu, table_hdu])
            return hdulist.writeto(filename, clobber=clobber, **kwargs)








