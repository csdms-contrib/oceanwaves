#!/usr/bin/env python

from oceanwaves import donelanwave, jonswapwave, yvwave


def main():
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


if __name__ == '__main__':
    main()
