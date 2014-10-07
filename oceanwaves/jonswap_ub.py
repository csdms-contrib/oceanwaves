import numpy as np
from oceanwaves import qkhfs

def jonswap_ub(Hsig,Tp,h):
    """Calculate ubr and Tbr from Hsig and Tp using JONSWAP spectrum

    Parameters
    ----------
    Hsig : float
        Significant wave height (m)
    Tp : float
        Peak period (s)
    h : float
        Water depth (m)

    Returns
    -------
    (ubr, Tbr) :
        ubr = representative bottom orbital velocity (m/s)
        Tbr = representative bottom wave period (s)
        
    """
    # Patricia Wiberg, UVa
    # Last modified 9 Mar 2007
    # Recoded in Python Sep 2014

    ubr = []
    Tbr = []
    dffp = 0.01
    ffp = np.arange(0.2,5+dffp,dffp)
    nt = len(Hsig)  
    fp = 1. / Tp

    for i in range(nt):
        m0 = Hsig[i] ** 2.0 / 16
        f = ffp * fp[i]
        df = dffp * fp[i]
        kh = qkhfs(2 * np.pi * f, h[i])
        gam = 3.3
        beta = 1.25
        xi = 0
        sig = list(ffp)
        sig = np.array([0.07 if j < 1 else 0.09 for j in sig])
        a = - ((ffp - 1) ** 2)
        b = 2.0 * sig ** 2.
        fdiv = a/b
        ee = np.exp(fdiv)
        t2 = gam ** ee
        t1 = -beta * (ffp ** -4)
        a = ffp ** xi
        b = ffp ** 5.
        fdiv2 = a/b
        sffpn = fdiv2 * np.exp(t1) * t2
        dnm = np.sum (sffpn * dffp)
        XJ = 1.0 / dnm
        a = m0 * fp[i] ** 4. * XJ
        b = f ** 5.
        fdiv3 = a/b
        sf = fdiv3 * np.exp(t1) * t2
        a = 2 * np.pi * f
        b = np.sinh(kh)
        fdiv4 = (a/b) ** 2
        su = fdiv4 * sf
        ubr = np.sqrt(2.0 * np.sum(su * df))
        fr = np.sum(f * su * df) / np.sum(su * df)
        fz = np.sqrt(np.sum(f ** 2 * su * df) / np.sum(su * df))
        Tbr = 1.0 / fr
        Tbz = 1.0 / fz
        ubr.append(ubr)
        Tbr.append(Tbr)
        
    return ubr, Tbr