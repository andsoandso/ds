# ds


A python module for 1d dynamical systems analysis.  It's for learning not for real work.

To install download the zip of this code to somewhere on your python path, or likewise `git clone https://github.com/andsoandso/ds.git` into your python path. As it is undergoing frequent updates, the git solution is the convenient one.

This code requires [numpy](http://www.numpy.org/) which can be (easily) installed for most if not all major platforms, see [here](http://www.scipy.org/scipylib/download.html) for the details.  

Requirements aside, this module is most easily used inside [ipython](http://ipython.org/)

---

This module is tested on OS 10.9, python 2.7.4, numpy 1.6.1, but should work on any recent unix-like system just fine; Let me know if that is not the case.  If you use windows I'd take patches but have no interesting in testing or supporting.  Windows, well, it just makes me sad.


# Intro

## In discrete time

Some quick examples of use.  Let's examine the function `2.5x(1-x)`, a logistic equation with a rate of 2.5.

Iterate our function

        >>> from ds.discrete import iterate
        >>> iterate(lambda x: (2.5*x)*(1-x), 0.6, 5)
        [0.6000000000000001, 0.6, 0.6000000000000001, 0.6, 0.6000000000000001]

        >>> iterate(lambda x: (2.5*x)*(1-x), 0, 5)
        [0.0, 0.0, 0.0, 0.0, 0.0]

        >>> # Between 0-1, 0.6 is stable fixed
        >>> iterate(lambda x: (2.5*x)*(1-x), 0.0001, 5)
        [0.000249975, 0.0006247812812484374, 0.0015609773239975975, 0.003896351684478907, 0.009702925320074413]
       
        >>> iterate(lambda x: (2.5*x)*(1-x), 0.7, 5)
        [0.5250000000000001, 0.6234375000000001, 0.586907958984375, 0.6061175166629255, 0.5968476816432343] 
    
        >>> # Above 1, or below 0 is unstable
        >>> iterate(lambda x: (2.5*x)*(1-x), 1.1, 5)
        [-0.27500000000000024, -0.8765625000000011, -4.1123107910156325, -52.55852708229812, -7037.393240357408]

        >>> iterate(lambda x: (2.5*x)*(1-x), -.1, 5)
        [-0.275, -0.8765624999999999, -4.1123107910156245, -52.558527082297935, -7037.393240357359] 

Find fixed points seeding from [-1, -.1, 0.01, .5, 1].

        >>> from ds.discrete import fixed_points
        >>> set([dis.fn.fixed_point(lambda x: (2.5*x)*(1-x), xs) for xs in [-1, -.1, 0.01, .5, 1]])
        set([0.0, 0.60000000000000009])

Are they stable? The neighborhood is 0.1.
        
        >>> from ds.discrete import is_stable
        >>> is_stable(lambda x: (2.5*x)*(1-x), 0, 0.1)
        (False, False)

        >>> is_stable(lambda x: (2.5*x)*(1-x), 0, 0.1)
        (True, True)

Print a pretty phase diagram.

        >>> from ds.vis import phase_diagram
        >>> phase_diagram((0,.6), ((False, False), (True, True)))

        ------<-----*----------->-----*-----<-----------------------
                    0                 0.6


## (More) fun iterating the logistic equation

Still using discrete maths. The general logistic function is `(r*x)*(1-x)`, implemented here using a lambda, with the parameters ordered so we can use partial.

[Partial](http://docs.python.org/2/library/functools.html#functools.partial) info.

The rate is 3.1, seed is 0.1

		>>> from functools import partial
		>>> iterate(partial(lambda r, x: (r*x)*(1-x), 3.1), 0.1, 5)
		[0.1, 0.2790000000000001, 0.6235929, 0.727646864715729, 0.6143484054538055]
		
Or to clean things up a bit

		>>> rate = 3.1; x0 = 0.1; iterate(partial(lambda r, x: (r*x)*(1-x), rate), x0, 5)
		[0.1, 0.2790000000000001, 0.6235929, 0.727646864715729, 0.6143484054538055]
				
		

## In continuous time

Unlike discrete, this branch uses the fn of interest's derivative.  The API and the signatures are similar to their discrete time equivalents.

I also introduce the `solver` submodule for solving differential equations.  To use solverfns, PFA (partial function application) should be used to normalize their signature to `solverfn(deriv, t, deriv_args=()))`.  

The raw signature for `ds.solver.euler` is `euler(deltat, deriv, t, deriv_args=())` but we do `partial(solver.euler, 2.0)` to normalize it to match the `solver()` signature. It sets the time-step to 2.  

Iterate the temperature equilibration differential eq `lambda x: 0.2*(20-x)` with a room temperature of 20.

        >>> from ds.continuous import iterate
		>>> from functools import partial
	    >>> iterate(lambda x: 0.2*(20-x), 1, 20, (), 
	    ....        partial(solver.euler, 2.0))
		[1, 8.600000000000001, 8.600000000000001, 13.16, 13.16, 15.896, 15.896, 17.5376, 17.5376, 18.522560000000002, 18.522560000000002, 19.113536, 19.113536, 19.4681216, 19.4681216, 19.68087296, 19.68087296, 19.808523775999998, 19.808523775999998, 19.8851142656, 19.8851142656]
				
We know a fixed point for `lambda x: 0.2*(20-x)` is 20.  Is it stable?
		
		>>> is_stable(lambda x: 0.2*(20-x), 20, 0.1, 
        ....		args=(), xtol=1e-4, maxiter=500, 
        ....		solverfn=partial(solver.euler, 2.0)) 
		(True, True)
		
Print a pretty phase diagram.

		>>> phase_diagram(xfix=(20.0,), xstable=((True, True), ), size=60, offset=12)


		------>-----*-----<-----------------------------------------
		            20.0

