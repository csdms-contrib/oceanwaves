import sys
from math import *
import numpy as np

from .qkhfs import qkhfs
from .ubspecpar import ubspecpar


def main(argv=None):
    if argv is None:
        argv = sys.argv

    #get input data
    print 'Enter file name with Hs, Tp, [h] [no header] (or return to enter values):'
    fnamestr = raw_input('')

    if fnamestr:
        wave_in_file = open(fnamestr, 'rU')
        hs = []
        tp = []
        h = []
        for line in wave_in_file:
            line = line.strip()
            columns = line.split()
            hsr = float(columns[0])
            tpr = float(columns[1])
            hs.append(hsr)
            tp.append(tpr)
            if len(columns) == 3:
                hr = float(columns[2])
                h.append(hr)
        if len(columns) == 2: 
            print 'Enter water depth (m):'
            hstr = raw_input('')
            h = float(eval(hstr))
        wave_in_file.close()
    else:
        print 'Enter signficant wave height (m):'
        hsstr = raw_input('')
        hs = [float(eval(hsstr))]
        print 'Enter peak wave period (s):'
        tpstr = raw_input('')
        tp = [float(eval(tpstr))]
        print 'Enter water depth (m):'
        hstr = raw_input('')
        h = [float(eval(hstr))]
    hs = np.array(hs)
    tp = np.array(tp)
    
    print 'Enter J for JONSWAP spectrum or D for Donnely spectrum (default):'
    spstr = raw_input('')
    if spstr: 
        specform = spstr
    else: 
        specform = 'D'

    #test data
    #tp=array([10,10,10,10], dtype=float)
    #hs=array([1,2,3,4], dtype=float)
    #h = 10.
    #specform = 'J'

    #run calculation
    [ubr,Tbr,itr] = ubspecpar(hs,tp,h,specform)

    #output results
    if fnamestr:
        ub_out_file = open('ubspecout.txt', 'w')
        for row in range(len(ubr)):
            print>>ub_out_file, ubr[row], Tbr[row]
        wave_in_file.close()
    else:
        print  'ubr = ',ubr[0], 'm/s'
        print 'Tbr = ', Tbr[0], 's'


if __name__ == "__main__":
    sys.exit(main())
