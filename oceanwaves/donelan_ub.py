import sys
import numpy as np
from oceanwaves import qkhfs

def donelan_ub(Hsig,Tp,h):
    """Calculate ubr and Tbr from Hsig and Tp using Donelan spectrum

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
        
    For Donelan spectrum, if number of iterations > 5, then
               calculation did not converge and will exit.
    """
    # Patricia Wiberg, UVa
    # Last modified 9 Mar 2007
    # Recoded in Python Sep 2014

    ubr = []
    Tbr = []
    g = 9.81
    dffp = 0.01
    ffp = np.arange(0.2,5+dffp,dffp)
    nt = len(Hsig)  
    fp = 1. / Tp

    for i in range(nt):
        m0 = Hsig[i] ** 2.0 / 16
        f = ffp * fp[i]
        df = dffp * fp[i]
        kh = qkhfs(2 * np.pi * f,h[i])
        itcnt = 0
        xi = 1
        tol = 0.001
        m0sfD = 0.0
        while abs((m0 - m0sfD) / m0 > tol):
            fpbar = (m0 * fp[i] ** 4. / (g ** 2 * 6.635e-06)) ** (1 / 0.7)
            gam = 1.7
            beta = 1.0
            if fpbar >= 0.05:
                if fpbar > 0.159: gam = 6.5 + 2.606 * np.log(fpbar)
                sig = 0.08 + 0.0013 * fpbar ** -3
            else:
                sig = 0.08 + 0.0013 * 0.05 ** -3
            a = - ((ffp - 1) ** 2.)
            b = 2.0 * sig ** 2.
            fdiv = a/b
            ee = np.exp(fdiv)
            t2 = gam ** ee
            t1 = -beta * (ffp ** -4.)
            a = ffp ** xi
            b = ffp ** 5.
            fdiv2 = a/b
            sffpn = fdiv2 * np.exp(t1) * t2
            dnm = np.sum (sffpn * dffp)
            XD = 1.0 / dnm
            a = m0 * fp[i] ** 4. * XD
            b = f ** 4. * fp[i]
            fdiv3 = a/b
            sf = fdiv3 * np.exp(t1) * t2
            m0sfD = np.sum(sf * df)
            itcnt = itcnt + 1
            if itcnt > 5:
                print 'Too many iterations'
                sys.exit()

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
