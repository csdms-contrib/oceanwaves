import sys
from math import *
import numpy as np


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

