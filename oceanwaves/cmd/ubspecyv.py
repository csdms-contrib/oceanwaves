#!/usr/bin/env python
import sys
from math import *
import numpy as np

from oceanwaves import yvwave, qkhfs, ubspecpar


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
    [hs,tp] = yvwave(u10,h,x)
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
