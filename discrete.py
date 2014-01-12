import numpy as np


def iterate(fn, x0, T, args=()):
    """Iterate fn T times, starting at x0.
    
    Parameters
    ----------
    fn : function
        A function whose first argument is x
    x0 : float
        Initial condition/seed 
    T : int
        The number of iterations to calculate   
    args : tuple, optional
        Extra arguments to `fn`.

    Return
    ------
    orbit - list
        The orbit/itinerary

    Examples
    --------
    >>> # The fixed points
    >>> iterate(lambda x: (2.5*x)*(1-x), 0.6, 20)
    >>> iterate(lambda x: (2.5*x)*(1-x), 0, 20)
    >>>
    >>> # Between 0-1, 0.6 is stable fixed
    >>> iterate(lambda x: (2.5*x)*(1-x), 0.0001, 20)
    >>> iterate(lambda x: (2.5*x)*(1-x), 0.7, 20)
    >>>
    >>> # Above 1, or below 0 is unstable
    >>> iterate(lambda x: (2.5*x)*(1-x), 1.1, 20)
    >>> iterate(lambda x: (2.5*x)*(1-x), -.1, 20)
    >>>
    >>> # Some assertions confirming the above facts
    >>> assert iterate(lambda x: (2.5*x)*(1-x), 0.6, 20)[0] == 0.6
    >>> assert iterate(lambda x: (2.5*x)*(1-x), 0.0001, 20)[-1] == 0.6
    >>> assert iterate(lambda x: (2.5*x)*(1-x), 0.99, 20)[0] == 0.6
    >>> assert iterate(lambda x: (2.5*x)*(1-x), 0, 20)[0] == 0
    >>> assert ds.fn.iterate(lambda x: (2.5*x)*(1-x), 1.1, 20) > 100000000000
    >>> assert ds.fn.iterate(lambda x: (2.5*x)*(1-x), -.1, 20) > 100000000000
    >>>
    >>> $ Confirm length of returned orbit
    >>> assert len(ds.fn.iterate(lambda x: (2.5*x)*(1-x), 0.6, 20)) == 20
    """
    
    # Initialize the orbit with x0 and x1
    orbit = [x0, fn(float(x0), *args)]
    
    # Iterate until t == T
    for t in range(1, int(T)):
        xt = fn(orbit[t - 1], *args)
        orbit.append(xt)
        
    return orbit


def fixed_point(fn, x0, args=(), xtol=1e-8, maxiter=500):
    """
    ----
    THIS CODE BORROWED FROM SCIPY:

    `http://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.fixed_point.html#scipy.optimize.fixed_point`

    I didn't want a full scipy dependency, instead this requires only
    numpy.

    NOTE: For consistency with other functions in thie module, 
        x0 must be float not array-like as in the scipy version.  
        Use a list comprehension to search many seeds.
    ----

    Find a fixed point of the function.

    Given a function of one or more variables and a starting point, find a
    fixed-point of the function: i.e. where ``func(x0) == x0``.

    Parameters
    ----------
    fn : function
        A function whose first argument is x
    x0 : float
        Fixed point of function.
    args : tuple, optional
        Extra arguments to `func`.
    xtol : float, optional
        Convergence tolerance, defaults to 1e-08.
    maxiter : int, optional
        Maximum number of iterations, defaults to 500.

    Notes
    -----
    Uses Steffensen's Method using Aitken's ``Del^2`` convergence acceleration.
    See Burden, Faires, "Numerical Analysis", 5th edition, pg. 80

    Examples
    --------
    >>> from scipy import optimize
    >>> def func(x, c1, c2):
    ....    return np.sqrt(c1/(x+c2))
    >>> c1 = 10
    >>> c2 = 2
    >>> optimize.fixed_point(func, 1.2, args=(c1,c2))
    1.6542491578567586
    """
    
    p0 = x0
    for iter in range(maxiter):
        p1 = fn(p0, *args)
        p2 = fn(p1, *args)
        d = p2 - 2.0 * p1 + p0
        if d == 0.0:
            return p2
        else:
            p = p0 - (p1 - p0)*(p1 - p0) / d
        if p0 == 0:
            relerr = p
        else:
            relerr = (p - p0)/p0
        if np.abs(relerr) < xtol:
            return p
        p0 = p

    msg = "Failed to converge after {0} iterations, value is {0}".format(
            maxiter, p)
    raise RuntimeError(msg)


def is_stable(fn, xfix, ep, args=(), xtol=1e-4, maxiter=500):
    """Is fn stable on the (plus, minus) sides of xfix?

    Parameters
    ----------
    fn : function
        A function whose first argument is x
    xfix : float
        The fixed point 
    ep : float
        The neighborhood to search for stability (ep > 0)
    args : tuple, optional
        Extra arguments to `fn`.
    xtol : float, optional
        Convergence tolerance, defaults to 1e-04.
    maxiter : int, optional
        Maximum number of iterations, defaults to 500.
    """

    if ep < 0:
        raise ValueError("ep must be positive")

    xps = []
    xms = []
    search_range = np.arange(0.01, ep, ep/10)
    for ep_x in search_range:
        xp = iterate(fn, xfix + ep_x, maxiter, *args)[-1]
        xm = iterate(fn, xfix - ep_x, maxiter, *args)[-1] ## Save last x
        xps.append(xp)
        xms.append(xm)

    p = np.all(np.abs(np.asarray(xps) - xfix) < xtol)
    m = np.all(np.abs(np.asarray(xms) - xfix) < xtol)

    return (p, m)


def phase_diagram(xfix=(), xstable=(), size=60, offset=12):
    """Display a (text) phase diagram

    Parameters
    ---------
    xfix : tuple
        The fixed points
    xstable : tuple (of tuples)
        Boolean stability tuples (from is_stable())
    size : int (size > 20)
        Width of the phase line in characters
    offset : int (offset > 2; size*.1 < offset < size*.25)
        Min/max offset of fixed points on the line
    """
    
    n = len(xfix)
    xfix = list(xfix)
    xstable = list(xstable)

    if n < 1:
        print("No fixed points")
        return None
    
    # Set some (hopefully) sane display limits
    if len(xfix) != len(xstable):
        raise ValueError("xfix and xstable must have same number of elements")
    if size < 20:
        raise ValueError("size must be greater than 20")
    if offset/float(size) < 0.1:
        raise ValueError("offset must be 10 percent of size")
    if offset/float(size) > 0.25:
        raise ValueError("offset must be less then 25 percent of size")
    if offset < 1:
        raise ValueError("offset must be > 1")        

    # Init a phase linees and its annotations
    line = ["-", ] * size
    annote = [" ", ] * size

    # Add fixed points
    steps = range(offset, size-offset, (size - 2*offset)/n)
    for s, xf in zip(steps, xfix):
        line[s] = "*"
        annote[s] = str(xf)

    # Add stablity arrows
    for i in range(len(xstable)):
        s = steps[i]
        xsb = xstable[i]

        downside = s - int(offset/2)
        upside = s + int(offset/2)

        # Plus side
        if xsb[0]:
            line[downside] = ">"
        elif not xsb[0]:
            # First step or last has no arrow?
            if (i == 0) or (not xstable[i-1][1]):
                line[downside] = "<"

        # Minus side
        if xsb[1]:
            line[upside] = "<"
        elif not xsb[1]:
            # Are we on the end?
            if i == len(xstable)-1:
                line[upside] = ">"

    print("\n")
    print(''.join(line))
    print(''.join(annote))
    
