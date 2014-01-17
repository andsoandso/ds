"""Bad ODE solvers.... Use scipy man, use scipy!"""


def euler(deltat, deriv, t, deriv_args=()):
    if deltat < 0:
        raise ValueError("deltat must be greater than 0")
    
    return (deriv(t, *deriv_args) * deltat, deltat)


def rk():
    raise NotImplementedError("TODO")


def adaptive():
    raise NotImplementedError("TODO")

