#! /usr/bin/env python

"""Wrappers for photutils' finder functions `find_peaks` and `daofind`
(technically `photutils/detection/findstars.daofind`).
I add on functionality to print output Table to a text file, similar to
IRAF's out files. The output file differs from IRAF, however, in that
it is column-based. It does include a parameter settings list in the
header comments.

Author:

    C.M. Gosmeyer, Feb. 2016

References:

About IRAF daofind:
http://stsdas.stsci.edu/cgi-bin/gethelp.cgi?daofind

About photutils daofind:
http://photutils.readthedocs.org/en/latest/_modules/photutils/detection/findstars.html#daofind

About photutils find_peaks:
http://photutils.readthedocs.org/en/latest/api/photutils.detection.find_peaks.html#photutils.detection.find_peaks

Example for use:
http://photutils.readthedocs.org/en/latest/photutils/getting_started.html

"""
from __future__ import print_function

import multiprocessing
import os
import sys
import photutils
import numpy as np
from astropy.io import fits
from astropy.stats import mad_std
from astropy.table import Column

from photutils_plus.phot_tools import read_in_photfile
from photutils_plus.phot_tools import write_out_photfile


#-------------------------------------------------------------------------------# 

def daofind(imagename, ext=1, outname='default', threshold=3., fwhm=1.5, \
            ratio=1.0, theta=0.0, sigma_radius=1.5, sharplo=0.2, \
            sharphi=1.0, roundlo=-1.0, roundhi=1.0, sky=0.0, \
            exclude_border=True, use_bkgrd=True):
    """ Extends `photutils.daofind`, so that outputs a coordinate file
    with my formatting from `photutils_plus.phot_tools.write_out_photfile`.
    (Later, if I'm clever, use inheritance.)

    Parameters
    ----------
    imagename : string
        Name of the FITS file.
    ext : int
        The extension of the FITS file to read.
    outname : string
        Name of the output coordinate file. If "default", becomes,
        "<imagename>.coo."
    threshold : int
        Threshold for search. Default of "3."
    fwhm : float
        The full width half max. Default of "1.5."
    ratio : float
        Default of "1.0."
    theta : float 
        Default of "0.0."
    sigma_radius : float
        Default of 1.5."
    sharplo : float
        Default of "0.2."
    sharphi : float
        Default of "1.0."
    roundlo : float
        Default of "-1.0."
    roundhi : float
        Default of "1.0."
    sky : float
        Default of "0.0."
    exclude_border : {True, Flase}
        When True, masks out a 10-pixel border on image from search.
    use_bkgrd : {True, False}
        Setting on will multiply the background by the threshold to
        obtain a new threshold value. 

    Returns
    -------
    coo_tab : astropy.Table
        Table containing the coordinates.

    """
    # Fetch function metadata.
    current_params = locals()
    func_name = sys._getframe().f_code.co_name #daofind.__name__

    # Read in FITS file.
    hdulist = fits.open(imagename)
    data = hdulist[ext].data

    # Mask 10 pixels from border. I find the `exclude_border` in 
    # photutils.daofind inadiquate. 
    if exclude_border:
        data[:,0:10] = -99999.0
        data[:,-10:] = -99999.0
        data[0:10,:] = -99999.0
        data[-10:,:] = -99999.0

    # Get the background.
    if use_bkgrd:
        bkg_sigma = mad_std(data) 
        threshold = threshold * bkg_sigma 

    coo_tab = photutils.daofind(data=data,
                                threshold=threshold,
                                fwhm=fwhm,
                                ratio=ratio,
                                theta=theta,
                                sigma_radius=sigma_radius,
                                sharplo=sharplo,
                                sharphi=sharphi,
                                roundlo=roundlo,
                                roundhi=roundhi,
                                sky=sky,
                                exclude_border=exclude_border)

    # ??Insert background value into coo_tab??

    # Create basename for out file.
    if outname == 'default':
        baseoutname = imagename + '.coo'
    else:
        baseoutname = outname

    # Check whether file already exists. If yes, append number.
    fileoutname = baseoutname
    i=1
    while os.path.isfile(fileoutname):
        fileoutname = baseoutname
        fileoutname += ('.' + str(i))
        i += 1

    write_out_photfile(fileoutname, coo_tab, current_params, func_name)

    return coo_tab


#-------------------------------------------------------------------------------# 

def find_peaks(imagename, ext=1, outname='default', threshold=40.0,
               box_size=24.0, border_width=5.0):
#, footprint=footprint, mask=mask, 
#               border_width=border_width, npeaks=npeaks, subpixel=subpixel,
#               error=error, wcs=wcs):
    """ Extends `photutils.find_peaks`.

    Parameters
    ----------
    imagename : string
        Name of the FITS file.
    ext : int
        The extension of the FITS file to read.
    outname : string
        Name of the output coordinate file. If "default", becomes,
        "<imagename>.coo."
    threshold : float
        Default of "40.0."
    box_size : float
        Default of "24.0."
    border_width : float
        Default of "5.0."
        
    Returns
    -------
    coo_tab : astropy.Table
        Table containing the coordinates.

    References
    ----------  
    http://photutils.readthedocs.org/en/latest/api/photutils.detection.find_peaks.html#photutils.detection.find_peaks
    
    """
    # Fetch function metadata.
    current_params = locals()
    func_name = sys._getframe().f_code.co_name 

    # Read in FITS file.    
    hdulist = fits.open(imagename)
    data = hdulist[ext].data
    hdulist.close()

    coo_tab = photutils.find_peaks(data=data, \
                                   threshold=threshold, \
                                   box_size=box_size, \
                                   border_width=border_width)
                                   #footprint=footprint, \
                                   #mask=mask, \
                                   #npeaks=npeaks, \
                                   #subpixel=subpixel, \
                                   #error=error, \
                                   #wcs=wcs)

    # Create basename for out file.
    if outname == 'default':
        baseoutname = imagename + '.coo'
    else:
        baseoutname = outname

    # Check whether file already exists. If yes, append number.
    fileoutname = baseoutname
    i=1
    while os.path.isfile(fileoutname):
        fileoutname = baseoutname
        fileoutname += ('.' + str(i))
        i += 1

    write_out_photfile(fileoutname, coo_tab, current_params, func_name)

    return coo_tab
     

#-------------------------------------------------------------------------------# 

def multiprocess_find(func, filenames):
    """ Wrapper for multiprocessing a finder function over a list of 
    FITS files.

    Parameters
    ----------
    func : find function
        Either :func:`daofind` or :func:`find_peaks`.
    filenames : list of strings
        The names of the FITS files to be multiprocssed with a 
        find function.
    """

    pool = multiprocessing.Pool(processes=SETTINGS['cores'])
    pool.map(func, filenames)
    pool.close()
    pool.join()

    print("Processes complete.")


#-------------------------------------------------------------------------------# 
#-------------------------------------------------------------------------------# 

def test_daofind():
    """
    """
    imagename = 'n6791_11924_55106_F606W_360_NOCI_F000_c0_1.fits' #'icqu04rgq_flc.fits' #'id1l02g2q_flt.clean.fits' 
    ext = 0
    daofind(imagename, ext=ext, threshold=1050., fwhm=2.5)

def test_find_peaks():
    """
    """
    imagename = 'id1l02g2q_flt.clean.fits'  #'icqu04rgq_flc.fits'
    ext = 0 
    find_peaks(imagename, ext=ext, threshold=40., box_size=24.)   


#-------------------------------------------------------------------------------# 
# The Main, for testing.
#-------------------------------------------------------------------------------# 

if __name__=='__main__':
    test_daofind()
    #data, header = read_in_photfile('id1l02g2q_flt.clean.fits.coo.1')
    #print data
    #print header




