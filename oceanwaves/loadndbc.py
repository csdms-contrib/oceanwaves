import numpy as np


def doy(month, day, year):
    """
    Calculation of day of year. If year is provided it will be tested for
    leap years.

    Parameters
    ----------
    month : numpy.ndarray or int32
        Month.
    day : numpy.ndarray or int32
        Day.
    year : numpy.ndarray or int32, optional
        Year.

    Retruns
    -------
    doy : numpy.ndarray or int32
        Day of year.
    """
    daysPast = np.array([0, 31, 60, 91, 121, 152, 182, 213, \
                         244, 274, 305, 335, 366])

    day_of_year = daysPast[month - 1] + day

    if year is not None:
        nonleap_years = np.invert(is_leap_year(year))
        day_of_year = day_of_year - nonleap_years + \
                      np.logical_and(day_of_year < 60, nonleap_years)

    return day_of_year


def is_leap_year(year):
    """
    Check if year is a leap year.

    Parameters
    ----------
    year : numpy.ndarray or int32

    Returns
    -------
    leap_year : numpy.ndarray or boolean
        True if year is a leap year.
    """
    return np.logical_or(np.logical_and(year % 4 == 0, year % 100 != 0),
                         year % 400 == 0)
    
def loadndbc():
    # LOADNBDC - Script to load NDBC text files
    # Works for new format with minutes in date

    # csherwood@usgs.gov
    # Last revised 18 Aug. 2006
    # Recoded in python by PL Wiberg Oct 2014

    #get input data
    print 'Enter file name NDBC spectral data:'
    fnamestr = raw_input('')

    ndbc_in_file = open(fnamestr, 'rU')
    hdr = ndbc_in_file.readline()
    hdr = hdr.strip()
    hcol = hdr.split()
    fstr = hcol[6:]
    nf = len(fstr)
    f=[]
    for i in range(nf):
       f.append(float(fstr[i]))
    yd = []        
    sr = []
    hs = []
    slong = []
    dvec = []
    rno = 1
    for line in ndbc_in_file:
        line = line.strip()
        columns = line.split()
        YYr = int(columns[0])
        MMr = int(columns[1])
        DDr = int(columns[2])
        HHr = int(columns[3])
        mnr = int(columns[4])
        sstr = columns[5:]
        for i in range(nf):
            sr.append(float(sstr[i]))
        hsr = 4*np.sqrt(np.sum(sr) / 100 )
        ydr =  float(doy(MMr, DDr, YYr)) + (float(HHr) +float(mnr) / 60.) / 24.
        dvecr = [YYr,MMr,DDr,HHr,mnr]
        yd.append(ydr)
        dvec.append(dvecr)
        hs.append(hsr)
        rno = rno + 1
    slong = np.array(sr)
    s = np.reshape(slong,(rno-1,nf))
    f = np.array(f)
    yd = np.array(yd)
    hs = np.array(hs)
        
    ndbc_in_file.close()

    # calculate signficant wave height (assumes uniform frequency bin size)
    
    return f,s,hs,yd

#[f,s,hs,yd] = loadndbc()
