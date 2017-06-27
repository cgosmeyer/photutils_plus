# photutils_plus
Extended functionality on the `astropy`-associated `photutils` photometry package.  

## Set Up

Clone the repo and run the setup.py.  You will need `astropy` and `photutils` installed.  

```
>>> python setup.py develop
```

I recommend "develop" rather than "install" because you'll probably want to tweak things.  (I'm not the Photometry Expert by a long shot.)  This was originally written in Python 2 and there may be some compatibility issues in Python 3 that I'm not aware of, but do open an issue if you find them!  

## Modules

### apphot
Can be run on the command line with limited amount of arguments (limited only by the free time of the developer) or be imported with full functionality.  See `apphot`'s doc strings for full details of the options. 

`apphot` wraps `photutils.CircularAperture` and `photutils.aperture_photometry` into a single function with parameters similar to those of IRAF's `DAOPHOT` task.  It outputs a photometry file with all the selected parameters saved in the header, including background estimates, similar to IRAF's `.mag` output files (but with a format easier to parse - see example at end of README).  In fact, the files are named as `.mag` by default to aid the transition for those used to their IRAF photometry. 

Still missing is the ability to re-centroid sources. To my knowledge, `photutils` has not developed their centroiding functions to make them able to handle images containing more than one source. 

Here are some command-line `apphot` examples, using an image named "id1l02g2q_flt.fits" and coordinate file named "id1l02g2q_flt.fits.coo.1" for aperture radii 5 and 10 pixels:

```
 >>> python apphot --im id1l02g2q_flt.fits --coo id1l02g2q_flt.fits.coo.1 --r 5.0 10.0
    
# This will generate two output *mag files in the current working directory named 
    
    id1l02g2q_flt.fits_5.0.mag
    id1l02g2q_flt.fits_10.0.mag
    
# If you set --sep 0, then there will be only one output *mag file named

    id1l02g2q_flt.fits.mag
    
# If the *mag files already exist, they will not be overwritten.
# Instead, the new *mag files will be appended with a number, i.e.,

    id1l02g2q_flt.fits_5.0.mag.1
    id1l02g2q_flt.fits_10.0.mag.1
```

### phot_tools
A module containing tools for the `apphot` and `sourcefind` modules.  Includes background estimating, masking, and file reading and writing. 

### plot_coords
Include in your photometry script to plot the coordinates found with `sourcefind` or another function onto a PNG of the FITS image. 

### sourcefind
Wraps `photutils.daofind` and more limitedly, `photutils.find_peaks`.  At this time, it cannot be run from the command-line (which would be easy to include if the developer took the initiative!).  

`sourcefind.daofind` outputs coordinate files with with a default naming scheme is like IRAF's: the image name appended with `.coo`.  

## A Scripted Workflow

For a single image, find the sources with `sourcefind.daofind`. 

```
from photutils_plus.sourcefind import daofind

imagename = 'icqu04rgq_flc.fits'
ext = 0
daofind(imagename, ext=ext, threshold=1050., fwhm=2.5)
```

This outputs a coordinate file named, by default, after the input image with ".coo" appended: "icqu04rgq_flc.fits.coo". 

Now, if you like, plot the coordinates as open circles onto the image.

```
from photutils_plus.plot_coords import plot_coords

plot_coords('icqu04rgq_flc.fits', 'icqu04rgq_flc.fits.coo', ext=1)
```

By default, the output PNG image will be named ""icqu04rgq_flc.fits.png".

Finally, perform photometry using `apphot.apphot`. By selecting "sep_out=False", we are confining all the different aperture radii (units of pixels) to a single file rather than a seperate file for each radius.  

```
from photutils_plus.apphot import apphot

apphot(imagename='icqu04rgq_flc.fits', ext=1, coofile='icqu04rgq_flc.fits.coo', ap_radii=[5, 10, 12], \
sep_out=False)
```

The output photometry file will be named, by default, as "icqu04rgq_flc.fits.mag". You can also have it return the magnitudes as an `astropy.table.Table` on top of the `.mag` file if you set `apphot` equal to a variable. This goes for `sourcefind.daofind` as well. 

## Example Output Photometry File

Here is an example of a `.mag` file from `apphot`, where "sep_out=False", so that all the aperture radii are printed in the same file for each source (3 shown in the example below). The header lists all the parameter settings; see the  `apphot` doc strings for more information on available parameters. 

```
# PYTHON                 2.7.5           version
# USER                   cgosmeyer               name
# HOST                   cgosmeyermac             computer
# DATE                   2016-08-13              yyyy-mm-dd
# TIME                   12:43:36                hh:mm:ss
# FUNCTION               apphot                  name
#
# ANNULUS                3.5             0
# AP_RADII               [3, 5, 10]            1
# BACKGLOBAL             False           2
# BACKMETHOD             mean            3
# COLNAMES               ['extr_xpix', 'extr_ypix']              4
# COOFILE                icqu04rgq_flc.coo                  5
# DANNULUS               10.0            6
# EFFECTIVE_GAIN                 None            7
# EXT            1               8
# IMAGENAME              numpy array             9
# MASK           None            10
# METHOD                 subpixel                11
# OUTNAME                icqu04rgq_flc           12
# PIXELWISE_ERROR                None            13
# SEP_OUT                False           14
# SUBPIXELS              5               15
# VERBOSE                False           16
#
ID      radius  aperture_sum    xcenter ycenter mean_local_bkgrd        tot_local_bkgrd   
1       3.0     304.709299699   177.559345818   2038.99871288   5.34541805349   151.138134785   
1       5.0     441.448582702   177.559345818   2038.99871288   2.35808650071   185.203680679    
1       10.0    823.806982231   177.559345818   2038.99871288   6.91177559481   2171.39834319   
2       3.0     253.913460445   213.889087553   2041.26538736   5.64445712521   159.593265342   
2       5.0     327.771958152   213.889087553   2041.26538736   12.0727819215   948.194074825    
2       10.0    11394.4246024   213.889087553   2041.26538736   2.13410972259   670.450342644   
3       3.0     619.62800065    28.2171662015   1995.3848863    5.61276484794   158.697187314   
3       5.0     795.264226607   28.2171662015   1995.3848863    3.13170050909   245.963182815    
3       10.0    1877.91108042   28.2171662015   1995.3848863    2.74348838443   861.892295375
...
```

