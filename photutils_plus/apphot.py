#! /usr/bin/env python

"""Wrapper for photutils' aperture photometry functions.
Meant to emulate some of the functionality of IRAF's DAOPHOT.

Author:

    C.M. Gosmeyer, Feb. 2016

Use:
    
    Can either be imported from or run from the command line.
    If from the command line, it has the following arguments

    --im - Name of the FITS file.

    --ext (optional) - Extension of FITS file to be read. Default of 1.

    --coo - Name of the coordinate file.

    --out (optional) - Name desired for output *mag file. By default, it
                       is <imagename>.mag<.n>, where 'n' increases from 1
                       if the *mag file already exists.

    --r (optional) - The aperture radii in units of pixels. Just 12.0 by 
                     default.

    --sep (optional) - 1 or 0. If on (1), will create a seperate out file 
                      for each aperture radius. If off (0) will add a row
                      for each aperture for each source. Off by default.
    
Example:

    >>> python apphot --im id1l02g2q_flt.fits --coo id1l02g2q_flt.fits.coo.1 --r 5.0 10.0

    This will generate two output *mag files in the cwd named 

    id1l02g2q_flt.fits_5.0.mag
    id1l02g2q_flt.fits_10.0.mag

    If you set --sep 0, then there will be only one output *mag file names

    id1l02g2q_flt.fits.mag

    If the *mag files already exist, they will not be overwritten.
    Instead, the new *mag files will be appended with a number, ie,

    id1l02g2q_flt.fits_5.0.mag.1
    id1l02g2q_flt.fits_10.0.mag.1   

References:

http://photutils.readthedocs.org/en/latest/photutils/aperture.html

http://www.astro.yale.edu/astr255/iraf/IRAF-daofind.html

Improvements:

    * Check that if *coo file was generated from IRAF, to then 
      subtract 1 pixel to x and y coords. (since IRAF 1-based)

"""

from __future__ import print_function

import os
import sys
import numpy as np
from astropy.io import fits
from astropy.io import ascii
from astropy import table
from astropy.table import Column
from astropy.table import Table
from photutils import aperture_photometry #as phot
from photutils import background
try: 
    from photutils import centroid_2dg
except:
    from photutils.morphology import centroid_2dg
from photutils import CircularAperture  # Aperture for pixel coords
from photutils import CircularAnnulus

from photutils_plus.phot_tools import meanclip_bkgrd
from photutils_plus.phot_tools import read_in_photfile
from photutils_plus.phot_tools import write_out_photfile

from photutils_plus.meanclip import meanclip


#-------------------------------------------------------------------------------# 

def fix_aperture_photometry(phot_tab, bad_phot_tab):
    """
    """
    phot_tab['aperture_sum'] = bad_phot_tab['aperture_sum']
    phot_tab['xcenter'] = (bad_phot_tab['xcenter'])[0][0]
    phot_tab['ycenter'] = (bad_phot_tab['ycenter'])[0][0]

    return phot_tab

#-------------------------------------------------------------------------------# 

def apphot(imagename, ext=1, coofile='', colnames='default', 
           outname='default',
           ap_radii=[12.0], sep_out=False, backmethod='mean',
           backglobal=False, method='subpixel', subpixels=5,
           centroid=False,
           pixelwise_error=None, effective_gain=None, mask=None,
           annulus=4., dannulus=2., verbose=False):
    """Performs aperture photometry for all given aperture radii (pixels).

    *mag files will be generated in cwd. They will not be overwriten. A 
    new file will be appended with '.n', where n is (1,2,3,...).


    ## Need a switch for doing either full-frame background or background
    ## via annulus?

    Parameters
    ----------
    imagename : string or array
        Either name of the FITS file or the array containing the data
        of that extension.
    ext : int
        Extension of FITS file to be read. Default of 1.
    coofile : string
        Name of the coordinate file.
    colnames : list of strings
        The column names of the x and y. If default, assumes they are
        'xcentroid' and 'ycentroid', which daofind returns.
        A non-default example is ['col3', 'col4'] (base of 1).
    outname : string
        Name desired for output *mag file. By default, it is
        <imagename>.mag<.n>, where 'n' increases from 1 if the *mag
        file already exists.
    ap_radii : list of floats or ints
        The apertures. Units of pixels.
    sep_out : {True, False}
        If on, will create a seperate out file for each aperture 
        radius. If off will add a row for each aperture for each
        source. Off by default.
    backmethod : string
        Method of calculating background. Currently 'mean' only
        supported option.
    backglobal : {True, False}
        Do either global background subtraction (where sources get
        masked and mean/median/mode taken of frame) or local (using
        the annulus and dannulus values). By default False - so 
        background is local by default.
    method : string
        Method to use for determining overlap between the aperture and 
        pixels. Options include ['center', 'subpixel', 'exact'], but not 
        all options are available for all types of apertures. More 
        precise methods will generally be slower.
    subpixels : int, optional
        If 'subpixel' is selected in method parameter, resample 
        pixels by this factor (in each dimension). That is, each 
        pixel is divided into subpixels ** 2 subpixels.
    centroid : {True, False}
        Turn on to centroid sources with ``phoutils.morphology.centroid_2dg``.
        Only recommended for single-source images for now.
    pixelwise_error : {True, False}
        For ``error`` and/or ``effective_gain`` arrays. If `True`,
        assume ``error`` and/or ``effective_gain`` vary significantly
        within an aperture: sum contribution from each pixel. If
        `False`, assume ``error`` and ``effective_gain`` do not vary
        significantly within an aperture. Use the single value of
        ``error`` and/or ``effective_gain`` at the center of each
        aperture as the value for the entire aperture.  Default is
        `True`.
    #error : 
    #effective_gain : 
    mask : array of bool
        Must be same dimensions as image.
    annulus : int or float
        The factor to multiply aperture radius by. For calculating 
        background.
    dannulus : int or float
        The depth of the annulus in pixels.
    verbose : {True, False}
        Set to True if want to print the tables to screen. False by default.

    Returns
    -------
    phot_tab or mega_phot_tab : `~astropy.table.Table`
        A table of the photometry with the following columns:

        * ``'aperture_sum'``: Sum of the values within the aperture.
        * ``'aperture_sum_err'``: Corresponding uncertainty in
          ``'aperture_sum'`` values.  Returned only if input ``error``
          is not `None`.
        * ``'xcenter'``, ``'ycenter'``: x and y pixel coordinates of the
          center of the apertures. Unit is pixel.
        * ``'xcenter_input'``, ``'ycenter_input'``: input x and y
          coordinates as they were given in the input ``positions``
          parameter.


    Outputs
    -------
    If 'sep_out' True, then for each radius in 'ap_radii',
    <imagename><_radius>.mag. If the file already exists, a
    number will be appended. An image named abcdefg.fits with 
    aperture radii of 5.0 and 12.0 will get *mag out files named, 
    for example, as

        abcdefg.fits_5.0.mag
        abcdefg.fits_5.0.mag.1
        abcdefg.fits_5.0.mag.2
        ...

    and

        abcdefg.fits_12.0.mag
        abcdefg.fits_12.0.mag.1
        abcdefg.fits_12.0.mag.2
        ...

    If 'sep_out' False, then all radii are in same file for each 
    source and the file is named as <imagename>.mag. Then naming
    scheme for an image name of abcdefg.fits is

        abcdefg.fits.mag
        abcdefg.fits.mag.1
        abcdefg.fits.mag.2
        ...

    In both cases of 'sep_out', <imagename> can be replaced by 
    'outname' if default is not set.

    The *mag files contain the five columns 'ID', 'radius', 
    'aperture_sum', 'xcenter', and 'ycenter'.
    Each row is a new ID. If 'sep_out' False then each
    ID will have a new row for each aperture.

    Notes
    -----
    Doesn't handle singletons well!  
    If only 1 star in coo file, it will NOT join tables (for some reason),
    so will need set sep_out=True. 
    Also, it will print the x and y coords to file within brackets
    (again, for some reason) and I don't know what to do about that.

    """
    # Fetch function metadata.
    current_params = locals()
    func_name = sys._getframe().f_code.co_name 

    # Read in FITS file.
    if '.fits' in imagename:
        hdulist = fits.open(imagename)
        data = hdulist[ext].data
        hdulist.close()
    else:
        data = imagename
        current_params['imagename'] = 'numpy array'

    # If a mask is given, set the parameter so shows nicely in outfile.
    if mask != None:
        current_params['mask'] = 'numpy array'

    # Read the coo file to get the positions.
    coo_tab = ascii.read(coofile)
    print(coo_tab)
    if len(coo_tab) == 0:
        print("Coordinate file is empty. Returning...")
        return Table()
    
    if colnames == 'default':
        xcoords = list(coo_tab['xcentroid'])
        ycoords = list(coo_tab['ycentroid'])
    else:
        xcoords = list(coo_tab[colnames[0]])
        ycoords = list(coo_tab[colnames[1]])
    positions = [xcoords, ycoords]
    
    # If centroiding coordinates set to True,
    if centroid:
        positions = centroid_2dg(data)

    # If finding a global background, only need do this once. 
    if backglobal:
        # This will add a 'global_bkgrd' column.
        # This does not quite work at this time if a source is near an edge.
        # Also, need have an exception if the imagename is a data array
        bkgrd, bkgrd_rms = meanclip_bkgrd(imagename, ext, backmethod, xcoords, ycoords)

    # Create mega-table if a single file for all apertures is selected.
    if not sep_out:
        mega_phot_tab = Table()

    for radius in ap_radii:
        radius = float(radius)
        # Get aperture objects using pixel coord positions.
        apertures = CircularAperture(positions, r=radius)

        # Do local background subtraction
        #annulus_apertures = CircularAnnulus(positions, r_in=radius*annulus, r_out=(radius*annulus)+dannulus)
        annulus_apertures = CircularAnnulus(positions, r_in=annulus, r_out=annulus+dannulus)
        # Get an error
        # this should be sky_sig?? how get that?

        # Then do photometry on the ap radius objects
        # the table has three columns, named 'aperture_sum', 'xcenter', 
        # and 'ycenter'.
        # The returned table is potentially 'bad' because aperture_photometry 
        # returns the x and y coords as numpy arrays if there is only one 
        # row of them. It only does this to coords, for reasons only it 
        # knows. This breaks stuff. We fix it below.
        bad_phot_tab = aperture_photometry(data, 
                                       apertures, 
                                       method=method,
                                       subpixels=subpixels, 
                                       pixelwise_error=pixelwise_error,
                                       #effective_gain=effective_gain,
                                       mask=mask)  #,error


        phot_tab = Table()

        if len(bad_phot_tab) == 1:
            # Change xcenter and ycenter so not in tuples.
            phot_tab = fix_aperture_photometry(phot_tab, bad_phot_tab)

        else:
            # Rows should be properly formatted; just set nice table
            # equal to the bad (now we know not-bad) table.
            phot_tab = bad_phot_tab

        # Find backgrounds, either local or global.
        if backglobal:
            # Add column to table. The global background was already 
            # calculated since it only needs be done once.
            phot_tab['global_bkgrd'] = bkgrd

        elif not backglobal:
            # Get a background table for annulus apertures.
            # This will add 'mean_local_bkgrd' and 'tot_local_bkgrd' columns,
            # where the first is the mean bkgrd per pixels in the annulus
            # and the total is the mean multiplied by the area of the source's aperture,
            # which should be used to subtract from the source's flux.
            bkgrd_phot_tab = Table()
            bad_bkgrd_phot_tab = aperture_photometry(data, 
                                                     annulus_apertures, 
                                                     method=method, 
                                                     subpixels=subpixels, 
                                                     pixelwise_error=pixelwise_error,
                                                     mask=mask)
            if len(bad_bkgrd_phot_tab) == 1:
                bkgrd_phot_tab = fix_aperture_photometry(bkgrd_phot_tab, 
                                                         bad_bkgrd_phot_tab)
            else:
                bkgrd_phot_tab = bad_bkgrd_phot_tab

            ### possibly make the below into a function ###
            if backmethod == 'mean':
                # Calculate mean local background
                bkgrd_mean = bkgrd_phot_tab['aperture_sum'] / annulus_apertures.area()
                bkgrd_sum = bkgrd_mean * apertures.area()
                phot_tab['mean_local_bkgrd'] = bkgrd_mean
                phot_tab['tot_local_bkgrd'] = bkgrd_sum

        # Add two more columns: one for ID and one for Aperture.
        col_id = Column(name='ID', data=np.arange(1, len(phot_tab['xcenter'])+1))
        col_ap = Column(name='radius', data=radius*np.ones(len(phot_tab['xcenter'])))
        phot_tab.add_column(col_id, 0)
        phot_tab.add_column(col_ap, 1)


        # Prepare to write the phot_table in a *mag file.

        # Create basename for out file.
        if outname == 'default':
            if sep_out:
                baseoutname = '{}_{}.mag'.format(coofile.split('.coo')[0], radius)
            else:
                baseoutname = '{}.mag'.format(coofile.split('.coo')[0])

        else:
            if sep_out:
                baseoutname = '{}_{}.mag'.format(outname, radius)
            else:
                baseoutname = '{}.mag'.format(outname)

        # Check whether file already exists. If yes, append number.
        fileoutname = baseoutname
        i=1
        while os.path.isfile(fileoutname):
            fileoutname = baseoutname
            fileoutname = '{}.{}'.format(fileoutname, i)
            i += 1
        
        if sep_out:
            # If creating a seperate *mag file for earch aperture.
            write_out_photfile(fileoutname, phot_tab, current_params, func_name)

        elif (not sep_out) and (len(mega_phot_tab) == 0):
            # If this is the first iteration, mega_phot_tab is empty and
            # can just be set equal to phot_tab.
            mega_phot_tab = phot_tab
            print("Created mega_phot_tab: {}".format(mega_phot_tab))

        else:
            # If this is a subsequent iteration, mega_phot_tab already exists and 
            # so just need add rows.
            print("Trying to join phot_tab to mega_phot_tab")
            if verbose:
                print(phot_tab)
                print(mega_phot_tab)
            mega_phot_tab = table.join(mega_phot_tab, phot_tab, join_type='outer')

    if not sep_out:
        # If writing all apertures to single file, write out the mega-table.
        write_out_photfile(fileoutname, mega_phot_tab, current_params, func_name)
        # Clear variable.
        del data
        
        return mega_phot_tab

    else:
        # Clear variable.
        del data

        return phot_tab



#-------------------------------------------------------------------------------#    
    
def parse_args():
    """Parses command line arguments.
    
    Returns
    -------
    args : object
        Containing the image and destination arguments.
            
    """

    im_help = 'Name of the FITS file.'
    ext_help = 'OPTIONAL. Extension of FITS file to be read. Default of 1.'
    coo_help = 'The coordinate file.'
    out_help = 'OPTIONAL. Custom name for out *mag file. Default is name ' + \
               'of the FITS file.'
    r_help = 'OPTIONAL. Aperture radii in pixels. Default of 12.0.'
    sep_help = 'OPTONAL. Set to 1 if you want seperate *mag files for ' + \
               'each aperture radius.'
        
    parser = argparse.ArgumentParser()
    parser.add_argument('--im', dest = 'imagename',
                        action = 'store', type = str, required = True,
                        help = im_help)
    parser.add_argument('--ext', dest = 'ext',
                        action = 'store', type = str, required = False,
                        help = ext_help, default=1)
    parser.add_argument('--coo', dest = 'coofile',
                        action = 'store', type = str, required = True,
                        help = coo_help)
    parser.add_argument('--out', dest = 'outname',
                        action = 'store', type = str, required = False,
                        help = out_help, default='default')
    parser.add_argument('--r', dest = 'ap_radii',
                        action = 'store', type = str, required = False,
                        help = r_help, default=['12.'])
    parser.add_argument('--sep', dest = 'sep_out',
                        action = 'store', type = str, required = False,
                        help = sep_help, default=0)   
    args = parser.parse_args()
     
    return args

# -----------------------------------------------------------------------------

def test_args(args):
    """ Ensures arguments are properly formatted.

    Inspired by `automated_scripts/cal_uvis_make_darks/cal_uvis_make_darks.py`
    by M. Bourque.

    Paramters
    ---------
    args : object
        Containing the image and destination arguments.

    Returns
    -------
    args : object
        Containing the image and destination arguments. Modified for
        proper formatting.

    """

    if args.sep_out == 0:
        args.sep_out = False
    else:
        args.sep_out = True

    return args

#-------------------------------------------------------------------------------# 
#-------------------------------------------------------------------------------# 

def test_wrapper():
    """
    """
    ap_radii = [1,2,3,4,5,6,7,8,9,10,12,14,16,18,20,24,28,32,36,40,45,50,55,60,65,70]
    #[5., 10., 12.]
    # IRAF coords
    iraf_coo = 'id1l02g2q_flt.clean_cnts.fits0.coo.1' 
    # Python coords
    python_coo = 'id1l02g2q_flt.clean.fits.coo'

    imagename = 'id1l02g2q_flt.clean.fits' #'n6791_11924_55106_F606W_360_NOCI_F000_c0_1.fits' #'icqu04rgq_flc.fits'  # 
    coofile =  python_coo #'n6791_11924_55106_F606W_360_NOCI_F000_c0_1.fits.coo' #'icqu04rgq_flc.fits.coo.2' #
    ext = 0
    apphot(imagename, ext=ext, coofile=coofile, ap_radii=ap_radii, sep_out=False, backglobal=True)

def test_parse():
    args = parse_args()
    args = test_args(args)

    imagename = args.imagename
    ext = args.ext
    coofile = args.coofile
    ap_radii = args.ap_radii
    sep_out = args.sep_out

    appot(imagename=imagename, ext=ext, coofile=coofile, ap_radii=ap_radii, sep_out=sep_out)

#-------------------------------------------------------------------------------# 
# The Main, for testing.
#-------------------------------------------------------------------------------# 

if __name__=='__main__':

    test_wrapper()


