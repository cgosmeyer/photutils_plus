# photutils_plus
Extended functionality on the `astropy`-associated `photutils` photometry package.  


## apphot
Can be run on the command line with limited amount of arguments (limited only by the free time of the developer) or be imported with full functionality.  

`apphot` wraps `photutils.CircularAperture` and `photutils.aperture_photometry` into a single function with parameters similar to those of IRAF's `DAOPHOT` task.  It outputs a photometry file with all the selected parameters saved in the header, including background estimates, similar to IRAF's `.mag` output files (but with a format easier to parse!).  In fact, the files are named as `.mag` by default to aid the transition for those used to their IRAF photometry. 

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

## 
