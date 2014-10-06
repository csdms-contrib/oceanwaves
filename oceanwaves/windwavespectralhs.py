import sys
from math import *
import numpy as np

def jonswapwave(u10,x,vflag):
    
#     JONSWAPWAVE(u10,x,vflag) returns Hsig,Tp [Tm,Tz]
#     Uses JONSWAP (Haseelmann et al 1973) wave spectrum 
#     to determine wave parameters Hsig and peak period
#
#     ALL UNITS MKS
#
#     INPUTS:
#     u10   = Windspeed at 10 m above surface (m/s)
#     x     = Wave fetch  (m)
#     vflag = c for constant (gam = 3.3) or v for veriable gamma spectrum
#                  
#     OUTPUTS:
#     Hsig  = Significant wave height (m)
#     Tp     = Peak period (s)
#     Tm     = Mean period (s)
#     Tz     = Zero-crossing period (s)
#
#     NOTES:
#     Hsig = 4*sqrt(int(S)), S= surface elevation spectrum
##
#     References:
#     Hasselman et al.

    g = 9.806
    dffp = 0.02
    ffp = np.arange(0.2,5+dffp,dffp)

    # Lower limits of u10,d,x
    minreal = 1.0e-4
    u10 = np.maximum(u10,minreal)
    x   = np.maximum(x, minreal)
    
    if x != 999:
        fp = 2.84 * g ** 0.7 * x ** -0.3 * u10 ** -0.4
        fpbar = fp * u10 / g
    else:
        fpbar = 0.13 
        fp = fpbar * g / u10
    f = ffp * fp    
    df = dffp * fp
    sig = list(ffp)
    sig = np.array([0.07 if j < 1 else 0.09 for j in sig])
    if vflag == 'c':
        alph = 0.033 * fpbar ** 0.667
        gam = 3.3
    else:
        alph = 0.033 * fpbar ** 0.86
        gam = 4.42 * fpbar ** 0.43
    beta = 1.25
    xi = 0
    eterm = -((ffp - 1) ** 2) / (2 * sig ** 2)
    ee = np.exp(eterm)
    t2 = gam ** ee
    t1 = -beta * (ffp ** -4)
    gp = g ** 2 / (2 * np.pi) ** 4
    sffp = alph * gp * (ffp ** xi / f **5) * np.exp(t1) * t2
    sffpn = (ffp ** xi / ffp **5) * np.exp(t1) * t2
    Xf = 1 / np.sum(sffpn * dffp)
    m0fp = np.sum(sffp * df)
    m1fp = np.sum(f * sffp * df)
    m2fp = np.sum(f **2. * sffp * df)
    Tzfp = np.sqrt(m0fp / m2fp)
    Tmfp = m0fp / m1fp

    # Calculate Hsig
    hs = 4 * np.sqrt(m0fp)
    
    Hsig = hs
    Tm = Tmfp
    Tz = Tzfp
    Tp = 1./fp
    
    return Hsig, Tp

def donelanwave(u10,x):
    
#     DONELANWAVE(u10,x) returns Hsig,Tp [Tm,Tz]
#     Uses Donelan (REF) wave spectrum 
#     to determine wave parameters Hsig and peak period
#
#     ALL UNITS MKS
#
#     INPUTS:
#     u10   = Windspeed at 10 m above surface (m/s)
#     x     = Wave fetch  (m)
#                  
#     OUTPUTS:
#     Hsig  = Significant wave height (m)
#     Tp     = Peak period (s)
#     Tm     = Mean period (s)
#     Tz     = Zero-crossing period (s)
#
#     NOTES:
#     Hsig = 4*sqrt(int(S)), S= surface elevation spectrum
##
#     References:
#    

    g = 9.806
    dffp = 0.02
    ffp = np.arange(0.2,5+dffp,dffp)

    # Lower limits of u10,d,x
    minreal = 1.0e-4
    u10 = np.maximum(u10,minreal)
    x   = np.maximum(x, minreal)
    
    if x != 999:
        fp = 2.84 * g ** 0.7 * x ** -0.3 * u10 ** -0.4
        fpbar = fp * u10 / g
    else:
        fpbar = 0.13 
        fp = fpbar * g / u10
    f = ffp * fp    
    df = dffp * fp
    gam = 1.7
    if fpbar >= 0.05:
        if fpbar > 0.159: gam = 6.5 + 2.606 * np.log(fpbar)
        sig = 0.08 + 0.0013 * fpbar ** -3
    else:
        sig = 0.08 + 0.0013 * 0.05 ** -3
    alph = 0.0165 * fpbar ** 0.55
    beta = 1.0
    xi = 1
    eterm = -((ffp - 1) ** 2) / (2 * sig ** 2)
    ee = np.exp(eterm)
    t2 = gam ** ee
    t1 = -beta * (ffp ** -4)
    gp = g ** 2 / (2 * np.pi) ** 4
    sffp = alph * gp * (ffp ** xi / f **5) * np.exp(t1) * t2
    sffpn = (ffp ** xi / ffp **5) * np.exp(t1) * t2
    Xf = 1 / np.sum(sffpn * dffp)
    m0fp = np.sum(sffp * df)
    m1fp = np.sum(f * sffp * df)
    m2fp = np.sum(f ** 2. * sffp * df)
    Tzfp = np.sqrt(m0fp / m2fp)
    Tmfp = m0fp / m1fp

    # Calculate Hsig
    hs = 4 * np.sqrt(m0fp)
    
    Hsig = hs
    Tm = Tmfp
    Tz = Tzfp
    Tp = 1./fp
    
    return Hsig, Tp


def yvwave(u10,d,x):
    
#     YVWAVE(u10,d,x) returns Hsig, Tpeak
#     Uses Young and Verhgen (1996) fetch-limited, finite-depth wave growth 
#     to determine wave parameters Hsig and peak period
#
#     ALL UNITS MKS
#
#     INPUTS:
#     u10   = Windspeed at 10 m above surface (m/s)
#     d     = Water Depth (m)
#     x     = Wave fetch  (m)
#                  
#     OUTPUTS:
#     Hsig  = Significant wave height (m)
#     Tp     = Peak period (s)
#
#     NOTES:
#     Hsig = 4*sqrt(E), E= variance(surface elevation)
#
#     u10 and d should be averaged along the fetch,e.g:
#         u10 = 1/x * integral([0->x],u10(xt)dxt) and
#         d   = 1/x * integral([0->x],d(xt)dxt)
#
# References:
#    Young, I.R. and L.A. Verhagen, 1996. The growth of fetch limited 
#     waves in water of finite depth. Part 1. Total energy and peak 
#     frequency. Coastal Engineering, 29: 47-78.

#     J.Paul Rinehimer 10 April 2007
#     Recoded in python by PL Wiberg, Sep 2014

    # Parameters
    minreal = 1.0e-4
    g       = 9.806
    n       = 1.74
    m       = -0.37
    invn = 1.0 / n
    invm = 1.0 / m
  
    # Lower limits of u10,d,x
    u10 = np.maximum(u10,minreal)
    d   = np.maximum(d, minreal)
    x   = np.maximum(x, minreal)
    
    # Calculate non-dimensional parameters
    chi   = g * x / u10 ** 2
    delta = g * d / u10 ** 2
  
    # Calculate a1,b1,a2,b2
    a1 = (0.292 * (delta ** 1.3)) ** invn
    b1 = (4.396E-5 * chi) ** invn
    a2 = (1.505 * (delta ** -0.375)) ** invm
    b2 = (16.391 * (chi ** -0.27)) ** invm
  
    ta1 = np.tanh(a1)
    ta2 = np.tanh(a2)
  
    # Calculate non-dimensional energy and frequency
    epsl = 0.00364 * (ta1 * np.tanh(b1/ta1)) ** n
    nu   = 0.133   * (ta2 * np.tanh(b2/ta2)) ** m
  
    # Calculate Hsig and period
    hs = 4 * np.sqrt(epsl * u10 ** 4 / g ** 2)
    t    = 1.0 / (g * nu / u10)
    
    Hsig = hs
    Tp = t
    
    return Hsig, Tp


print 'Enter wind speed 10m above surface (m/s): '
u10str = raw_input('')
u10 = float(eval(u10str))
print 'Enter J (JONSWAP), D (Donelan) or Y (Young&Verhagen) spectrum [default = D]:'
sflag = raw_input('')
if sflag: sflag = sflag.upper()
if sflag == '': sflag = 'D'
if sflag == 'J' or sflag == 'D':
    print 'Enter fetch length in km [default = fully developed sea]: '
    xstr = raw_input('')
    if xstr:
        x = float(eval(xstr))
        x = x * 1000.
    else:
        x = 999
if sflag == 'D': [Hsig,Tp] = donelanwave(u10,x)
elif sflag == 'J':
    print 'Enter c (gam=const) or m (variable gamma) JONSWAP spectrum [default = m]:'
    vflag = raw_input('')
    if vflag == '': vflag = 'm'
    [Hsig,Tp] = jonswapwave(u10,x,vflag)
elif sflag == 'Y': 
    print 'Enter water depth [default = 2m]:'
    dstr = raw_input('')
    if dstr:
        d = float(eval(dstr))
    else:
        d = 2.
    print 'Enter fetch length in km [default = 10 km]: '
    xstr = raw_input('')
    if xstr:
        x = float(eval(xstr))
        x = x * 1000.
    else:
        x = 10000.
    [Hsig,Tp] = yvwave(u10,d,x)
else:
    print 'Rerun using a valid spectrum identifier (J,D or Y)'
    sys.exit()

    
print 'Hsig = ',Hsig,' m'
print 'Tpeak = ',Tp,' s; fpeak = ',1./Tp,' s-1'
