ds
==

A python module for dynamical systems analysis.  It's for learning not for real work.

To install download the zip of this code to somewhere on your python path, or likewise `git clone https://github.com/andsoandso/ds.git` into your python path.

This code requires [numpy](http://www.numpy.org/) which can be (easily) installed for most if not all major platforms, see [here](http://www.scipy.org/scipylib/download.html) for the details.

Intro
=====

Some quick examples of use.  Let's examine the function `2.5x(1-x)`.

Iterate our function

        >>> from ds.fn import iterate
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

Find the fixed pointsi seeding from [-1, -.1, 0.01, .5, 1].

        >>> from ds.fn import fixed_points
        >>> set(dis.fn.fixed_point(lambda x: (2.5*x)*(1-x), [-1, -.1, 0.01, .5, 1]))
        set([0.0, 0.60000000000000009])

Are they stable (within a neighborhood of 0.1)?
        
        >>> from ds.fn import is_stable
        >>> is_stable(lambda x: (2.5*x)*(1-x), 0, 0.1)
        (False, False)

        >>> is_stable(lambda x: (2.5*x)*(1-x), 0, 0.1)
        (True, True)

Print a pretty phase diagrami to the console.

        >>> ds.fn.phase_diagram((0,.6), ((False, False), (True, True)))

        ------------*----------->-----*-----<-----------------------
                    0                 0.6
