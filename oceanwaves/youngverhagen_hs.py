import numpy as np


def youngverhagen_hs(u10,d,x):
    """Significant wave height and peak period using Young and Verhgen spectrum.
    
    Uses Young and Verhgen (1996) fetch-limited, finite-depth wave growth 
    to determine wave parameters Hsig and peak period
    
    Parameters
    ----------
    u10 : float
        Windspeed at 10 m above surface (m/s)
    d : float
        Water Depth (m)
    x : float
        Wave fetch (m)

    Returns
    -------
    (Hsig, Tp)
        Significant wave height (m) and peak period (s)

    Notes
    -----
    Hsig = 4*sqrt(E), E = variance (surface elevation)

    u10 and d should be averaged along the fetch,e.g:
        u10 = 1/x * integral([0->x],u10(xt)dxt) and
        d   = 1/x * integral([0->x],d(xt)dxt)

    References
    ----------
    Young, I.R. and L.A. Verhagen, 1996. The growth of fetch limited 
        waves in water of finite depth. Part 1. Total energy and peak 
        frequency. Coastal Engineering, 29: 47-78.
    """
    # J.Paul Rinehimer 10 April 2007
    # Recoded in python by PL Wiberg, Sep 2014

    # Parameters
    minreal = 1.0e-4
    g       = 9.806
    n       = 1.74
    m       = -0.37
    invn = 1.0 / n
    invm = 1.0 / m
  
    # Lower limits of u10,d,x
    u10 = np.maximum(u10,minreal)
    d   = np.maximum(d, minreal)
    x   = np.maximum(x, minreal)
    
    # Calculate non-dimensional parameters
    chi   = g * x / u10 ** 2
    delta = g * d / u10 ** 2
  
    # Calculate a1,b1,a2,b2
    a1 = (0.292 * (delta ** 1.3)) ** invn
    b1 = (4.396E-5 * chi) ** invn
    a2 = (1.505 * (delta ** -0.375)) ** invm
    b2 = (16.391 * (chi ** -0.27)) ** invm
  
    ta1 = np.tanh(a1)
    ta2 = np.tanh(a2)
  
    # Calculate non-dimensional energy and frequency
    epsl = 0.00364 * (ta1 * np.tanh(b1/ta1)) ** n
    nu   = 0.133   * (ta2 * np.tanh(b2/ta2)) ** m
  
    # Calculate Hsig and period
    hs = 4 * np.sqrt(epsl * u10 ** 4 / g ** 2)
    t    = 1.0 / (g * nu / u10)
    
    Hsig = hs
    Tp = t
    
    return Hsig, Tp
