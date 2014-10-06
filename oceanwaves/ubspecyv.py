import sys
from math import *
import numpy as np
import yvwave

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
  


def ubspecpar(hs,tp,h,specform):
    # UBSPECPAR - Caclulate ubr and Tbr from hs and tp using parametric spectrum
    # [ubr,Tbr,iter] = ubspecform(hs,tp,h,specform)
    #
    # Input:
    #   hs - Significant wave height (m)
    #   tp - Peak period (s)
    #   h  - Water depth (m)
    #   specform - spectral formulation to use
    #       specform='D' for Donelan spectrum (default)
    #       specform='J' for JONSWAP spectrum
    # Returns:
    #   ubr = representative bottom orbital velocity (m/s)
    #   Tbr = representative bottom wave period (s)
    #   iter = number of iterations if Donelan spectrum is chosen
    #       If iter>5 the calculation did not converge
    #   iter=0 for the JONSWAP spectrum
    #
    # Patricia Wiberg, UVa
    # Last modified 9 Mar 2007
    # Recoded in Python Sep 2014

    ubr = []
    Tbr = []
    itr = []
    g = 9.81
    dffp = 0.01
    ffp = np.arange(0.2,5+dffp,dffp)
    nt = np.size(hs)  
    fp = np.divide(1.0,tp)

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
 
def main(argv=None):
    if argv is None:
        argv = sys.argv

    #get input data
    print 'Enter file name with u10, fetch, depth [no header] (or return to enter values):'
    fnamestr = raw_input('')

    if fnamestr:
        wind_in_file = open(fnamestr, 'rU')
        u10 = []
        x = []
        h = []
        for line in wind_in_file:
            line = line.strip()
            columns = line.split()
            u10r = float(columns[0])
            xr = float(columns[1])
            hr = float(columns[2])
            u10.append(u10r)
            x.append(xr)
            h.append(hr)
        wind_in_file.close()
    else:
        print 'Enter wind at 10m (m/s):'
        u10str = raw_input('')
        u10 = [float(eval(u10str))]
        print 'Enter fetch (m):'
        xstr = raw_input('')
        x = [float(eval(xstr))]
        print 'Enter water depth (m):'
        hstr = raw_input('')
        h = [float(eval(hstr))]
    u10 = np.array(u10)
    x = np.array(x)
    h = np.array(h)
    print h, type(h)
    
    print 'Enter J for JONSWAP spectrum or D for Donnely spectrum (default):'
    spstr = raw_input('')
    if spstr: 
        specform = spstr
    else: 
        specform = 'D'

    #run Young and Verhagen model to get hs, tp
    [hs,tp] = yvwave.yvwave(u10,h,x)
    if type(hs) == list: 
        hs = np.array(hs)
    if hs.ndim == 2:
        hs = hs[0,:]
    if type(tp) == list: 
        tp = np.array(tp)
    if tp.ndim == 2:
        tp = tp[0,:]
    #print 'Hs = ',hs
    #print 'Tp = ',tp
    
    #run ubspec calculation
    [ubr,Tbr,itr] = ubspecpar(hs,tp,h,specform)

    #output results
    if fnamestr:
        ub_out_file = open('ubspecyvout.txt', 'w')
        for row in range(len(ubr)):
            print>>ub_out_file, ubr[row], Tbr[row]
    else:
        print  'ubr = ',ubr[0], 'm/s'
        print 'Tbr = ', Tbr[0], 's'
       
if __name__ == "__main__":
    sys.exit(main())

