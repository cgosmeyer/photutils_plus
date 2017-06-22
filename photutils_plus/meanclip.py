from __future__ import print_function

import numpy as np


def meanclip(indata, clipsig=3.0, maxiter=5, converge_num=0.02, verbose=1, return_array=0, return_median=0):
    """Computes an iteratively sigma-clipped mean on a data set. 
    
    Clipping is done about median, but mean is returned by default 
    (use return_median to return the clipped median value).
    
    Parameters:
        indata : array_like
            Input data.
        
        clipsig : float
            Number of sigma at which to clip.
        
        maxiter : int
            Ceiling on number of clipping iterations.
        
        converge_num : float
            If the proportion of rejected pixels is less than
            this fraction, the iterations stop.
        
        verbose : {0, 1}
            Print messages to screen?
        
        return_array : {0, 1}
            Return the final array indices that were used to compute
            statistics.
            
        return_median : {0, 1}
            Return the median if 1. Return the mean if 0. 
        
    Returns:
        val : float
            The N-sigma clipped mean or median.
        sigma : float
            The standard deviation of remaining pixels.
        
    Outputs:
        Prints to screen mean or median statistics.      
        
    History:
        * 21/10/1998 Written by RSH, RITSS
        * 20/01/1999 Added SUBS, fixed misplaced paren on float call, improved doc. RSH
        * 24/11/2009 Converted to Python. PLL.
        * 08/01/2013 Added option to return the array indices of non-clipped pixels. DMH
        * Added option to return median of the clipped array. DMH
    
    Notes:
        This is based on MYMEANCLIP routine from ACS library.    
    
    Examples:
    
    >>> mean, sigma = meanclip(indata)
    """
    # Flatten array
    skpix = indata.reshape( indata.size, )
    
    # initialize array to store indices of pixels used to compute stats
    arrind = np.arange(0,skpix.size)
    
    ct = indata.size
    iter = 0; c1 = 1.0 ; c2 = 0.0
    
    while (c1 >= c2) and (iter < maxiter):
        lastct = ct
        medval = np.median(skpix)
        sig = np.std(skpix)
        wsm = np.where( abs(skpix-medval) < clipsig*sig )
        ct = len(wsm[0])
        if ct > 0:
            skpix = skpix[wsm]
            arrind = arrind[wsm]
        
        c1 = abs(ct - lastct)
        c2 = converge_num * lastct
        iter += 1
    
    if return_median:
        val = np.median(skpix)
        val_type = 'median'
    else:
        val  = np.mean( skpix )
        val_type = 'mean'
    sigma = np.std( skpix )
    
    if verbose:
        if return_median:
            prf = 'MEDIANCLIP:'
            print('{} {}.1f-sigma clipped median'.format(prf, clipsig))
            print('{} Median computed in {} iterations'.format(prf, iter))
            print('P{} Median = {}, sigma = {}'.format(prf, val, sigma))
        
        else:
            prf = 'MEANCLIP:'
            print('{} {}.1f-sigma clipped mean'.format(prf, clipsig))
            print('{} Mean computed in {} iterations'.format(prf, iter))
            print('{} Mean = {}, sigma = {}'.format(prf, val, sigma))
    
    if return_array:
        return np.copy(arrind)
    else:
        return val, sigma