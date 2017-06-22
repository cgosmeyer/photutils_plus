
"""
Module for Python photomtery assisting functions.

Author: 

    C.M. Gosmeyer, March 2016
"""


import getpass
import numpy as np
import os
import platform 
import shutil
import socket
import time

from astropy.io import ascii
from astropy.io import fits
from collections import OrderedDict
from photutils_plus.meanclip import meanclip


#-------------------------------------------------------------------------------# 

def append_data_to_header(filename):
    """Appends data to header.
    
    Parameters
    ----------
    filename : string
        Name of the file.

    """
    with open(filename, "a") as header_file:
        data_file = open(filename + '.temp', "r")
        data_list = data_file.readlines()
        data_file.close()
        header_file.write('\n')
        for line in data_list:
            if '#' not in line:
                header_file.write(line)

    # Delete the temp file.
    os.remove(filename+'.temp')

#-------------------------------------------------------------------------------#

def circular_mask(arr_shape, r, x_offset=0, y_offset=0):
    """Generates circular mask for 2D image.
        
    Parameters
    ----------
    arr_shape : tuple of int
        Shape of the array to use the mask.
    r : int
        Radius of the mask in pixels.
    x_offset, y_offset : int or float, optional
        Mask offset relative to image center.

    Returns
    -------
    Numpy indices of the mask, rounded to nearest integer.
          
    References
    ----------
    http://mail.scipy.org/pipermail/numpy-discussion/2011-January/054470.html
    """
    assert len(arr_shape) == 2, 'Image is not 2-D'
    
    ny, nx = arr_shape
    assert nx > 1 and ny > 1, 'Image is too small'
    
    assert isinstance(r, (int, long)) and r > 0, 'Radius must be int > 0'
    
    xcen = np.round(0.5 * nx - 0.5 + x_offset).astype('int')
    ycen = np.round(0.5 * ny - 0.5 + y_offset).astype('int')
    
    x1, x2 = xcen - r, xcen + r
    y1, y2 = ycen - r, ycen + r

    assert y1 >= 0 and y2 < ny and x1 >= 0 and x2 < nx, \
          'Mask falls outside image bounds'
    
    y, x = np.ogrid[-r:r, -r:r]
    i = np.where(x**2 + y**2 <= r**2)
    
    a = np.zeros(arr_shape).astype('bool')
    a[y1:y2, x1:x2][i] = True
    
    return np.where(a)


#-------------------------------------------------------------------------------# 

def multiple_source_mask(arr_shape, r, xcoords, ycoords):
    """ Returns an array of True and Falses.
    The True are are where the (square!) masks are applied.
    Includes edge-detection!

    Parameters
    ----------
    arr_shape : array
        The width and height of the image array, [w,h].
    r : float or integer
        The radius of the mask.
    xcoords : array
        Array of the x positions.
    ycoords : array
        Array of the y positions.

    Returns
    -------
    Array of booleans.

    """
    ny, nx = arr_shape

    bool_array = np.zeros(arr_shape).astype('bool')

    for xc, yc in zip(xcoords, ycoords):
        x1, x2 = xc-r, xc+r
        y1, y2 = yc-r, xc+r

        # Make sure mask falls within the image. If it doesn't, trim.
        if y1 <= 0:
            y1 = 1
        elif y2 > ny:
            y2 = ny-1
        elif x1 <= 0:
            x1 = 1
        elif x2 > nx: 
            x2 = nx-1
        
        bool_array[y1:y2, x1:x2] = True

    return np.where(bool_array)

#-------------------------------------------------------------------------------# 

def meanclip_bkgrd(imagename, ext=0, backmethod='mean', xcoords=[], ycoords=[],
    naxis1='', naxis2='',detector='', maskrad='', multiple_sources=False):
    """ After masking out sources and image border, calculates
    the mean clipped background.

    Based off of D. Hammer's script from 
    `detectors.uvis_contam.ptscr_photom_uvis.run_daophot_uvis`
    """
    # Read in FITS file.
    if '.fits' in imagename:
         # Create temporary image for bckgrd measurement that masks sources 
        # out to 80 pixels (assign a very low number).
        tmp_imagename = imagename+'.back.fits'
        shutil.copy(imagename, tmp_imagename)
        hdulist = fits.open(tmp_imagename, mode='update')
        naxis1 = hdulist[0].header['naxis1']
        naxis2 = hdulist[0].header['naxis2']
        detector = hdulist[0].header['detector']
        maskim = hdulist[ext].data
    else:
        maskim = imagename


    if detector == 'IR' and maskrad == '':
        maskrad = 30
    elif detector == 'UVIS' and maskrad == '':
        maskrad = 80
    
    if multiple_sources:
        maskim[multiple_source_mask(maskim.shape, maskrad, xcoords, ycoords)] = -99999.0
    else:
        for xc, yc in zip(xcoords, ycoords):
            maskim[circular_mask(maskim.shape, maskrad, x_offset=(xc-naxis1/2.0), \
               y_offset=(yc-naxis2/2.0))] = -99999.0

    # Also mask out sources with zero effective exposure 
    # [WE ELIMINATE PIXELS WITHIN 20 OF IMAGE BORDER].
    maskim[:,0:20] = -99999.0
    maskim[:,-20:] = -99999.0
    maskim[0:20,:] = -99999.0
    maskim[-20:,:] = -99999.0

    # Generate initial guess for lower/upper limits (use 10 sigma).
    fmaskim = np.ndarray.flatten(maskim)
    llim = -100
    ulim = 10000.0
    init_median,init_rms = meanclip(fmaskim[(fmaskim > llim) & \
                                    (fmaskim < ulim)], maxiter=7, \
                                    return_median=1)
    llim = init_median - 10.0*init_rms
    ulim = init_median + 10.0*init_rms

    # Measure background and rms.

    if backmethod.lower() == 'mean':
        back,backrms=meanclip(fmaskim[(fmaskim > llim) & \
                              (fmaskim < ulim)], maxiter=7)
    elif backmethod.lower() == 'median':
        back,backrms = meanclip(fmaskim[(fmaskim > llim) & \
                                (fmaskim < ulim)], maxiter=7, 
                                return_median=1)
    elif backmethod.lower() == 'mode':
        backmean,backrms = meanclip(fmaskim[(fmaskim > llim) & \
                                    (fmaskim < ulim)], maxiter=7)
        nbins = np.ceil(80.0/(0.1*backrms))
        cc,bb,pp = pylab.hist(fmaskim[(fmaskim > llim) & \
                              (fmaskim < ulim)], log=True, bins=nbins, \
                              range=(-40.0,40.0))
        back = bb[cc.argmax()] + (bb.max()-bb.min())/(2.0*(len(bb)-1))
    else:
        raise Exception('Background statistical method {} is not' + \
                        ' covered in our case list.'.format(backmethod))

    return back, backrms


#-------------------------------------------------------------------------------# 

def read_in_photfile(filename, isheader=True):
    """ Reads in my custom photometry file, created in 
    :func:`write_out_photfile`.

    Parameters
    ----------
    filename : string
        Name of the file.
    isheader : {True, False}
        Set to True if the file contains a header.
    """
    data = ascii.read(filename)

    if isheader:
        header = ascii.read(data.meta['comments'], delimiter='\t',
                        format='no_header', names=['key', 'val', 'units'])
    else:
        header=None

    return data, header


#-------------------------------------------------------------------------------# 

def parse_phot_header(header):
    """ Parses my photometry file header.

    Parameters
    ----------
    header : dictionary
        The header of photometry file.
    """
    keys = header['key']
    vals = header['val']
    
    parsed_header = {}

    for key, val in zip(keys, vals):
        parsed_header[key.lower()] = val


    return parsed_header


#-------------------------------------------------------------------------------# 

def write_out_photfile(fileoutname, tab, current_params={}, func_name='', 
    header=False, verbose=False):
    """ Writes output photometry file with all parameters listed at 
    top in comments.

    Parameters
    ----------
    fileoutname : string
        Name of the output photometry file.
    tab : astropy.Table
        Table of the photometry outputs.
    current_params : dictionary
        The parameters used in the photometry function: ap_radii, 
        centroid, etc.
    func_name : string
        Name of the photometry function used. 
    header : {True, False}
        Set to True if header (and, presumably, the file) already exists.
    verbose : {True, False}
        Set to True for verbose mode.

    References
    ----------
    Fetching a function's parameters.
    http://stackoverflow.com/questions/582056/getting-list-of-parameter-names-inside-python-function

    Astropy ascii.read. See especially 'Comments and metadata' section.
    http://docs.astropy.org/en/stable/io/ascii/read.html

    Astropy ascii.write can't append. drrr
    https://github.com/astropy/astropy/issues/3684

    """
    # If already have a header, just insert it. If not, create it.

    if not header:
        # Order the current_params dictionary
        current_params = OrderedDict(sorted(current_params.items(), key=lambda t: t[0]))

        # Create logistic values to be inserted at top of file's header.
        python_version = platform.python_version()  # or sys.version
        user = getpass.getuser()
        host = socket.gethostname()
        date = time.strftime("%Y-%m-%d")
        hms = time.strftime("%H:%M:%S")

        # Define lists of logistics.
        logistic_keys = ['PYTHON', 'USER  ', 'HOST  ', 'DATE  ', 'TIME  ', 'FUNCTION']
        logistic_values = [python_version, user, host, date, hms, func_name]
        logistic_units = ['version', 'name', 'computer', 'yyyy-mm-dd', 'hh:mm:ss', 'name']


        # Open the file. 
        open_file = open(fileoutname, 'w')
        if verbose:
            print("Opening {} to be written...".format(open_file))

        # First write logistics.
        for key, val, unit in zip(logistic_keys, logistic_values, logistic_units):
            open_file.write("# {} \t\t {} \t\t {}\n".format(key, val, unit) )
        open_file.write("#\n")

        # Second write the finder function's parameters. 
        units = np.arange(len(current_params.keys()))
        for key, val, unit in zip(current_params.keys(), current_params.values(), units):
            #print("# {} \t\t {} \t\t {}\n".format(key.upper(), val, unit))
            open_file.write("# {} \t\t {} \t\t {}\n".format(key.upper(), val, unit) )

        open_file.write("#\n")

    else:
        # If input a header, write it to file.
        # Open the file. 
        open_file = open(fileoutname, 'w')
        for i in range(len(header[header.colnames[0]])):
            row_str = '# '
            for col in header[i]:
                row_str += str(col) + "\t"
            row_str += "\n"
            open_file.write(row_str)


    # Write out the column names of the table.
    colname_str = ''
    if verbose:
        print("printing colnames...")
    for colname in tab.colnames:
        colname_str += colname + "\t"
        if verbose:
            print(colname_str)
    colname_str += "\n"

    open_file.write(colname_str)

    # Finally write out the table, row by row.
    if verbose:
        print('printing table row by row...')
    for i in range(len(tab[tab.colnames[0]])):
        row_str = ''
        for col in tab[i]:
            col = str(col).split(' pix')[0]  # one day unit started appearing after all the coords. drr astropy tables drr
            row_str += str(col) + "\t"
        row_str += "\n"
        if verbose:
            print(row_str)
        open_file.write(row_str)

    # Close the file.
    open_file.close()
    

