def binarize(x, treshold = 1):
    """
    Accepts list of numbers x and shrinks numbers on interval [0,1].
    All numbers below treshold are mapped to 0.
    """
    maximal = max(x)
    return [0 if el<treshold else el/maximal for el in x]