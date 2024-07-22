import numpy as np

def signOP(x):
    """
    Custom sign function that returns -1 for negative and +1 for non-negative values.

    Parameters:
    x (numeric or numpy.ndarray): Input number or array.

    Returns:
    int or numpy.ndarray: The sign of the input as defined by the function.
    """
    s = 2 * (x >= 0) - 1
    return s

