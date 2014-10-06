import sys
from math import *
import numpy as np
import loadndbc as lx

def qkhfs(w,h):
    #QKHFS - Quick iterative calculation of kh in dispersion relationship
    # kh = qkhf( w, h )
    #
    # Input:
    #  w Angular wave frequency = 2*pi/T where T = wave period [1/s]
    #  h Water depth [m]
    # Returns: kh = wavenumber * depth [ ]
    # 
    # Either w or h can be a vector, but not both.
    # Hard-wired for MKS units.
    # Orbital velocities from kh are accurate to 3e-12 !
    #
    # RL Soulsby (2006) "Simplified calculation of wave orbital velocities"
    # HR Wallingford Report TR 155, February 2006
    # Eqns. 12a - 14

    # csherwood@usgs.gov
    # Sept 10, 2006
    # Recoded in Python by P. Wiberg Sep 2014

    g = 9.80665
    
    if type(w) == list:
        x = [item ** 2 * h / g for item in w]
    elif type(h) == list:
        x = [w ** 2 * item / g for item in h]
    else:
        x = w ** 2 * h / g
        
    lx = len(x)
    kh = [float(i) for i in x]
    kh = [sqrt(i) if i < 1 else i for i in x ]
    for i in range(lx):
        y = kh[i]
        t = tanh(y)
        y = y - ((y * t - x[i]) / (t + y * (1. - t ** 2)))
        t = tanh(y)
        y = y - ((y * t - x[i]) / (t + y * (1. - t ** 2)))
        t = tanh(y)
        y = y - ((y * t - x[i]) / (t + y * (1. - t ** 2)))
        kh[i] = y
    return kh
  


def ubspecdat(h,s,f,df):
# UBSPECDAT - Calculate ubr and Tbr from measured spectra
# ubspecdat(h,s,f,df) returns ubr,Tbr
#
# Input:
#   h = water depth (m) - scalar or col. vector with length
#   s(nf) or s(nt,nf) = array of spectral densities normalized so that
#                       Hs = 4*sum(s,2)*df
#   f       = row vector with central frequencies (Hz)
#   df      = (optional) scalar or row vector with freq. bandwidths (Hz)
# Returns:
#   ubr = representative bottom orbital velocity (m/s)
#   Tbr = representative bottom wave period [ (s)
#   The alternative bottom period, Tbz, is also calculated (see text).

# Chris Sherwood, USGS
# Last revised September 8, 2006
# Recoded in Python by PL Wiberg, Oct 2014

    print h
    print f
    print s    
    xx = s.shape
    print xx
    nf = xx[1]
    nt = xx[0]
    w = 2 * np.pi * f
    # Determine kh using Soulsby (2006) method 
    kh = qkhfs(w,h)
    w = np.tile(w,(nt,1))
    fm = np.tile(f,(nt,1))
    kh = np.tile(kh,(nt,1))
    h = h * np.ones((nt,nf))

    print df
    if df:
        if np.len(df)==1:
            df = df * np.ones(nt,nf)
        elif np.len(df) == nf:
            df = np.tile(df,(nt,1))
    else:
        df = np.diff(f)
        print type(df)
        print type(df[-1])
        df = np.hstack((df,df[-1]))
        df = np.tile(df,(nt,1))

    
    Su = ((w ** 2) / (np.sinh(kh) ** 2)) * s
    ubr = np.sqrt(2 * np.sum((Su * df)))
    fr = np.sum((Su * fm * df)) / (np.sum((Su * df)))
    Tbr = 1. / fr
    fz = np.sqrt(np.sum((Su * fm ** 2 * df))) / (np.sum((Su * df)))
    Tbz = 1. / fz
    
    return ubr, Tbr

[f,s,hs,yd] = lx.loadndbc()

h = 10
[ubr,Tbr] = ubspecdat(h,s,f,[])



