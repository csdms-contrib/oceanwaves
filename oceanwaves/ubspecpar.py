import numpy as np


def ubspecpar(hs,tp,h,specform):
    """Calculate ubr and Tbr from hs and tp using parametric spectrum

    Parameters
    ----------
    hs : float
        Significant wave height (m)
    tp : float
        Peak period (s)
    h : float
        Water depth (m)
    specform - {'D', 'J'}
        spectral formulation to use
            specform='D' for Donelan spectrum (default)
            specform='J' for JONSWAP spectrum

    Returns
    -------
    (ubr, Tbr, iter) :
        ubr = representative bottom orbital velocity (m/s)
        Tbr = representative bottom wave period (s)
        iter = number of iterations if Donelan spectrum is chosen
               If iter > 5 the calculation did not converge
               iter=0 for the JONSWAP spectrum
    """
    # Patricia Wiberg, UVa
    # Last modified 9 Mar 2007
    # Recoded in Python Sep 2014

    ubr = []
    Tbr = []
    itr = []
    g = 9.81
    dffp = 0.01
    ffp = np.arange(0.2,5+dffp,dffp)
    nt = len(hs)  
    fp = 1. / tp

    for i in range(nt):
        m0 = hs[i] ** 2.0 / 16
        f = ffp * fp[i]
        df = dffp * fp[i]
        kh = qkhfs(2 * pi * f,h[i])
        itcnt = 0
        if specform == 'D':
            xi = 1
            tol = 0.001
            m0sfD = 0.0
            while abs((m0 - m0sfD) / m0 > tol):
                fpbar = (m0 * fp[i] ** 4. / (g ** 2 * 6.635e-06)) ** (1 / 0.7)
                gam = 1.7
                if fpbar >= 0.05:
                    if fpbar > 0.159: gam = 6.5 + 2.606 * log(fpbar)
                    sig = 0.08 + 0.0013 * fpbar ** -3
                else:
                    sig = 0.08 + 0.0013 * 0.05 ** -3
                a = - ((ffp - 1) ** 2.)
                b = 2.0 * sig ** 2.
                fdiv = a/b
                ee = np.exp(fdiv)
                t2 = gam ** ee
                t1 = -(ffp ** -4.)
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

            a = 2 * pi * f
            b = np.sinh(kh)
            fdiv4 = (a/b) ** 2
            su = fdiv4 * sf
            ubrD = np.sqrt(2.0 * np.sum(su * df))
            frD = np.sum(f * su * df) / np.sum(su * df)
            fzD = np.sqrt(np.sum(f ** 2 * su * df) / np.sum(su * df))
            TbrD = 1.0 / frD
            TbzD = 1.0 / fzD
            ubr.append(ubrD)
            Tbr.append(TbrD)
            itr.append(itcnt)
        
        if specform == 'J':
            gam = 3.3
            xi = 0
            sig = list(ffp)
            sig = np.array([0.07 if j < 1 else 0.09 for j in sig])
            a = - ((ffp - 1) ** 2)
            b = 2.0 * sig ** 2.
            fdiv = a/b
            ee = np.exp(fdiv)
            t2 = gam ** ee
            t1 = -1.25 * (ffp ** -4)
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
            a = 2 * pi * f
            b = np.sinh(kh)
            fdiv4 = (a/b) ** 2
            su = fdiv4 * sf
            ubrJ = np.sqrt(2.0 * np.sum(su * df))
            frJ = np.sum(f * su * df) / np.sum(su * df)
            fzJ = np.sqrt(np.sum(f ** 2 * su * df) / np.sum(su * df))
            TbrJ = 1.0 / frJ
            TbzJ = 1.0 / fzJ
            ubr.append(ubrJ)
            Tbr.append(TbrJ)
            itr.append(itcnt)
        
    return ubr, Tbr, itr
