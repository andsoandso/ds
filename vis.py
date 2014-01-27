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
        raise ValueError("xfix and xstable must match")
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
    

def fixed_point(xfix=(), size=60, offset=12):
    """Display a (text) fixed point diagram

    Parameters
    ---------
    xfix : tuple
        The fixed points
    size : int (size > 20)
        Width of the phase line in characters
    offset : int (offset > 2; size*.1 < offset < size*.25)
        Min/max offset of fixed points on the line
    """
    
    n = len(xfix)
    xfix = list(xfix)

    if n < 1:
        print("No fixed points")
        return None
    
    # Set some (hopefully) sane display limits
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

    print("\n")
    print(''.join(line))
    print(''.join(annote))