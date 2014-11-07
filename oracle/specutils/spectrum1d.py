# coding: utf-8

""" Class for dealing with 1D spectra. """

from __future__ import division, print_function

__author__ = "Andy Casey <arc@ast.cam.ac.uk>"

__all__ = ["Spectrum1D", "Spectrum"]

# Standard library
import logging
import os

# Third-party
import numpy as np
import pyfits

logger = logging.getLogger("oracle")

class Spectrum1D(object):
    """
    This is a temporary class holder for a Spectrum1D object until the
    :class:`astropy.specutils.Spectrum1D` module has advanced sufficiently to replace it.
    """
    
    def __init__(self, disp, flux, variance=None, headers=None):
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
            # Assumed to be Poisson
            logger.warn("No variance given. Assuming noise to be Poisson-distributed.")
            self.variance = self.flux.copy()
        else:
            self.variance = variance.copy()

        self.ivariance = 1.0/self.variance
        if headers is not None:
            self.headers = headers
        else:
            self.headers = {}
        return None


    def copy(self):
        """ Creates a copy of the object """
        return self.__class__(self.disp.copy(), self.flux.copy(),
            variance=self.variance.copy(), headers=headers.copy())
    

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
            image = pyfits.open(filename, **kwargs)
            
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
                #li = a.headers['LTM1_1'] * np.arange(a.headers['NAXIS1']) + a.headers['LTV1']
                #a.headers['CRVAL1'] + a.headers['CD1_1'] * (li - a.headers['CRPIX1'])

                if np.all([key in header.keys() for key in ('CDELT1', 'NAXIS1', 'CRVAL1')]):
                    disp = header['CRVAL1'] + np.arange(header['NAXIS1']) * header['CDELT1']
            
                if "LTV1" in header.keys():
                    disp -= header['LTV1'] * header['CDELT1']

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
            try:
                disp, flux, variance = np.loadtxt(filename, unpack=True, **kwargs)
            except:
                disp, flux = np.loadtxt(filename, unpack=True, **kwargs)
            
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
            disp = pyfits.Column(name='disp', format='1D', array=self.disp)
            flux = pyfits.Column(name='flux', format='1D', array=self.flux)
            var = pyfits.Column(name='variance', format='1D', array=self.variance)
            table_hdu = pyfits.new_table([disp, flux, var])

            # Create Primary HDU
            hdu = pyfits.PrimaryHDU()

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
            hdulist = pyfits.HDUList([hdu, table_hdu])
            return hdulist.writeto(filename, clobber=clobber, **kwargs)



