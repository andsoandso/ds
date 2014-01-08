import numpy as np


def iterate(fn, x0, T, args=()):
    """Iterate fn T times, starting at x0.
    
    Parameters
    ----------

    fn : function
        A function that takes a single argument
    x0 : float
        Initial condition/seed 
    T : int
        The number of iterations to return   
    args : tuple, optional
        Extra arguments to `fn`.

    Return
    ------

    orbit - list
        The orbit/itinerary

    Examples/Tests
    --------------

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
    
    # Initialize the orbit with x0
    orbit = [fn(float(x0), *args), ]
    
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
    numpy.  I'd like to remove even that eventually.
    ----


    Find a fixed point of the function.

    Given a function of one or more variables and a starting point, find a
    fixed-point of the function: i.e. where ``func(x0) == x0``.

    Parameters
    ----------
    fn : function
        Function to evaluate.
    x0 : array_like
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
    >>> c1 = np.array([10,12.])
    >>> c2 = np.array([3, 5.])
    >>> optimize.fixed_point(func, [1.2, 1.3], args=(c1,c2))
    array([ 1.4920333 ,  1.37228132])
    """

    if not np.isscalar(x0):
        x0 = np.asarray(x0)
        p0 = x0
        for iter in range(maxiter):
            p1 = fn(p0, *args)
            p2 = fn(p1, *args)
            d = p2 - 2.0 * p1 + p0
            p = np.where(d == 0, p2, p0 - (p1 - p0)*(p1 - p0) / d)
            relerr = np.where(p0 == 0, p, (p-p0)/p0)
            if np.all(np.abs(relerr) < xtol):
                return p
            p0 = p
    else:
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


def is_stable(fn, xfix, dx, args=(), xtol=1e-4, maxiter=500):
    """Is fn stable on the (plus, minus) sides of xfix?

    Parameters
    ----------

    fn : function
        A function that takes a single argument
    xfix : float
        The fixed point 
    dx : float
        Amount to preturb xfix by
    args : tuple, optional
        Extra arguments to `fn`.
    xtol : float, optional
        Convergence tolerance, defaults to 1e-08.
    maxiter : int, optional
        Maximum number of iterations, defaults to 500.
    """

    p = False
    m = False

    xp = iterate(fn, xfix + dx, maxiter, *args)[-1]
    xm = iterate(fn, xfix - dx, maxiter, *args)[-1]
        ## Save only the last element

    if abs(xp - xfix) < xtol:
        p = True
    if abs(xm - xfix) < xtol:
        m = True

    return (p, m)


def phase_diagram(fn, minval, maxval):
    """Display as a (text) phase diagram"""
    
    raise NotImplementedError("TODO")

    fixed = find_fixed(fn)
    for fix in fixed:
        if is_stable(fn, fix):
            pass
        elif is_unstable(fn, fix):
            pass    
    pass

