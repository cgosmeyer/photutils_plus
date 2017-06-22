#! /usr/bin/env python

""" Plots the coordinates from a file onto a FITS image.

Author:

    C.M. Gosmeyer, Mar. 2016

"""

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np
from astropy.io import fits
from photutils_plus.phot_tools import read_in_photfile

#-------------------------------------------------------------------------------# 

def plot_coords(imagename, coofile, ext=1, xcol='default', ycol='default', 
    outname='default'):
    """ Overplot coordinates onto FITS image and saves as PNG.

    Parameters
    ----------
    imagename : string
        Name of the FITS file.
    coofile : string
        Name of the coordinate file.
    ext : int
        The extension of the FITS file to read.
    xcol : string
        The name of the x-coordinate column in the 'coofile.'
        If "default," then assumes column name is "xcentroid."
    ycol : string
        The name of the y-coordinate column in the 'coofile.'
        If "default," then assumes column name is "ycentroid."
    outname : string
        Name of the output PNG image.
    """
    if type(imagename) == str:
        # # Get image data from FITS file.
        hdulist = fits.open(imagename)
        imdata = hdulist[ext].data
        hdulist.close()
    else:
        # Or if the input is already image data,
        imdata = imagename

    # Get the x and y pixel coords.
    coodata, cooheader = read_in_photfile(coofile, isheader=False)

    if xcol == 'default': 
        xcoords = list(coodata['xcentroid'])
    else: 
        xcoords = list(coodata[xcol])

    if ycol == 'default':
        ycoords = list(coodata['ycentroid'])
    else:
        ycoords = list(coodata[ycol])

    # To adjust the contrast of image.
    #plt.hist(data.ravel(), bins=256,  fc='k', ec='k')
    #plt.show()

    shape = np.shape(imdata)
    plt.figure(figsize=(shape[0]/100., shape[1]/100.))

    plt.imshow(imdata, cmap='gray', origin='lower', clim=(-20, 40)) #, norm=LogNorm())
    plt.scatter(xcoords, ycoords, s=10, facecolors='none', edgecolors='r')
    #plt.show()
    
    if outname == 'default':
        diagplotname = '{}.png'.format(imagename)
    else:
        diagplotname = '{}.png'.format(outname)
    
    plt.savefig(diagplotname, bbox_inches='tight')

    return diagplotname
    

#-------------------------------------------------------------------------------# 
# The Main, for testing.
#-------------------------------------------------------------------------------# 
   
if __name__=="__main__":
    # Examples
    #plot_coords('icqu04rgq_flc.fits', 'icqu04rgq_flc.fits.coo.2', ext=1)
    #plot_coords('n6791_11924_55106_F606W_360_NOCI_F000_c0_1.fits', 'n6791_11924_55106_F606W_360_NOCI_F000_c0_1.fits.coo', ext=0)
    #plot_coords('n6791_11924_55106_F606W_360_NOCI_F000_c0_1.fits', '../orig_mag_coo/n6791_11924_55106_F606W_360_NOCI_F000_c0_1.fits.coo.1', ext=0)
