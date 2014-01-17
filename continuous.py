import numpy as np


def iterate(deriv, x0, T, args=(), solverfn=None):
    """Find the orbit (of length T) for x0 given its derivative.
    
    Parameters
    ----------
    derivative : function
        A derivative function whose first argument is t
    x0 : float
        Initial condition/seed 
    T : int
        The number of iterations to calculate   
    args : tuple, optional
        Extra arguments to `deriv`.
    solverfn : function
        A solver function with a signature matching: 
            `solverfn(deriv, x, t, args)`
    """
    
    if solverfn == None:
        raise ValueError("solverfn was not defined")
    
    # Initialize the orbit with x0
    t = 0
    orbit = [x0, ]
    
    # Iterate until i == T
    for i in range(0, int(T)):
        dx, dt = solverfn(deriv, orbit[i - 1], *args)
        xt = orbit[i - 1] + dx
        orbit.append(xt)
    
    return orbit


def is_stable(deriv, xfix, ep, args=(), xtol=1e-4, maxiter=500, solverfn=None):
    """Is deriv stable on the (plus, minus) sides of xfix?

    Parameters
    ----------
    deriv : function
        A derivative function whose first argument is t
    xfix : float
        The fixed point 
    ep : float
        The neighborhood to search for stability (ep > 0)
    args : tuple, optional
        Extra arguments to `deriv`.
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
        xp = iterate(deriv, xfix + ep_x, maxiter, args, solverfn)[-1]
        xm = iterate(deriv, xfix - ep_x, maxiter, args, solverfn)[-1] 
            ## Save last x
        xps.append(xp)
        xms.append(xm)
        
    p = np.all(np.abs(np.asarray(xps) - xfix) < xtol)
    m = np.all(np.abs(np.asarray(xms) - xfix) < xtol)

    return (p, m)


if __name__ == '__main__':
    from functools import partial
    from ds import solver

    print("Testing `iterate()`")
    print("\n")    
    
    orbits10 = iterate(lambda x: 0.2*(20-x), 1, 20, (), 
            partial(solver.euler, 2.0))
    print("Using `lambda x: 0.2*(20-x)` at T = 10, x0 = 0, and a dt of 2.0")
    print("the first 10 orbits were:")
    print(orbits10)
    print("\n")
    
    orbits100 = iterate(lambda x: 0.2*(20-x), 1, 200, (), 
            partial(solver.euler, 0.2))
    print("Using `lambda x: 0.2*(20-x)` at T = 100, x0 = 0, and a dt of .2")
    print("the first 100 orbits were:")
    print(orbits100)
    print("\n")
    
    print("Both should approach 20. If they don't warnings will follow.")
    assert abs(orbits10[-1] - 20) < 1, "orbits10 is off"
    assert abs(orbits100[-1] - 20) < 1, "orbits100 is off"
    print("Done.")
    
    print("Testing `is_stable()`")
    print('Using `lambda x: 0.2*(20-x)` an ep = 0.1,'
     ' at dt = 10, xfix = 20, and a dt of 2.0.')
    stab = is_stable(lambda x: 0.2*(20-x), 20, 0.1, 
            args=(), xtol=1e-4, maxiter=500, 
            solverfn=partial(solver.euler, 2.0))
    print(stab)
    assert stab == (True, True), "`is_stable()` is off."
    print("\n")
    
    print("The matching phase line.")
    from ds.vis import phase_diagram
    phase_diagram(xfix=20.0, xstable=stab, size=60, offset=12)
    