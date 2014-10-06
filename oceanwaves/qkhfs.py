
def qkhfs(w,h):
    """Quick iterative calculation of kh in dispersion relationship

    Parameters
    ----------
    w : float
        Angular wave frequency = 2 * pi / T where T = wave period [1/s]
    h : float
        Water depth [m]

    Returns
    -------
    kh : float
        wavenumber * depth [ ]

    Notes
    -----
    Either w or h can be a vector, but not both.
    Hard-wired for MKS units.
    Orbital velocities from kh are accurate to 3e-12 !

    References
    ----------
    RL Soulsby (2006) "Simplified calculation of wave orbital velocities"
    HR Wallingford Report TR 155, February 2006
    Eqns. 12a - 14
    """
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
