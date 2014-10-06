import sys
from math import *
import numpy as np

def jonswapwave(u10,x,vflag):
    """Significant wave height and peak period using jonswap method.

    Uses JONSWAP (Haseelmann et al 1973) wave spectrum 
    to determine wave parameters Hsig and peak period

    ALL UNITS MKS

    Parameters
    ----------
    u10 : float
        Windspeed at 10 m above surface (m/s)
    x : float
        Wave fetch  (m)
    vflag :{'c', 'v'}
        *c* for constant (gam = 3.3) or *v* for veriable gamma spectrum

    Returns
    -------
    (Hsig, Tp, Tm, Tz)
        Hsig = Significant wave height (m)
        Tp = Peak period (s)
        Tm = Mean period (s)
        Tz = Zero-crossing period (s)

    Notes
    -----
    Hsig = 4*sqrt(int(S)), S= surface elevation spectrum

    References
    ----------
    Hasselman et al.
    """
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
    """Significant wave height and peak period using Donelan method.

    Uses Donelan (REF) wave spectrum to determine wave parameters Hsig and
    peak period

    ALL UNITS MKS

    Parameters
    ----------
    u10 : float
        Windspeed at 10 m above surface (m/s)
    x : float
        Wave fetch  (m)

    Returns
    -------
    (Hsig, Tp, Tm, Tz)
        Hsig  = Significant wave height (m)
        Tp     = Peak period (s)
        Tm     = Mean period (s)
        Tz     = Zero-crossing period (s)

    Notes
    -----
    Hsig = 4*sqrt(int(S)), S = surface elevation spectrum
    """

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
