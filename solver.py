"""Bad ODE solvers.... Use scipy man, use scipy!"""


def euler(deltat, deriv, x, t, deriv_args=()):
    if deltat < 0:
        raise ValueError("deltat must be greater than 0")
    
    return x + deriv(t, *deriv_args) * deltat


def rk():
    pass


def adaptive():
    pass

