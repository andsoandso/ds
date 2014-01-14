

def iterate(deriv, x0, T, args=(), solverfn=None):
    """Find the orbit (of length T) for x0 given its derivative.
    
    Parameters
    ----------
    derivative : function
        A function whose first argument is t
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
        dx, dt = solverfn(deriv, x, t, *args)
        xt = orbit[i - 1] + dx
        orbit.append(xt)
        t += dt
    
    return orbit


