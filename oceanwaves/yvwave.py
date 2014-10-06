import sys
from math import *
import numpy as np

def yvwave(u10,d,x):
    
# [Hsig,t] = YVWAVE(u10,d,x)
#     Uses Young and Verhgen (1996) fetch-limited, finite-depth wave growth 
#     to determine wave parameters Hsig and peak period
#
#     ALL UNITS MKS
#
#     INPUTS:
#     u10   = Windspeed at 10 m above surface (m/s)
#     d     = Water Depth (m)
#     x     = Wave fetch  (m)
#                  
#     OUTPUTS:
#     Hsig  = Significant wave height (m)
#     t     = Peak period (s)
#
#     NOTES:
#     Hsig = 4*sqrt(E), E= variance(surface elevation)
#
#     u10 and d should be averaged along the fetch,e.g:
#         u10 = 1/x * integral([0->x],u10(xt)dxt) and
#         d   = 1/x * integral([0->x],d(xt)dxt)
#
# References:
#    Young, I.R. and L.A. Verhagen, 1996. The growth of fetch limited 
#     waves in water of finite depth. Part 1. Total energy and peak 
#     frequency. Coastal Engineering, 29: 47-78.

#     J.Paul Rinehimer 10 April 2007

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
    Hsig = []
    Tp = []
    
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
    
    Hsig.append(hs)
    Tp.append(t)
    
    return Hsig, Tp
