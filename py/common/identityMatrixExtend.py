import numpy as np

def identityMatrixExtend(Z, cyclicValue):
    """
    Transforms base graph to H by cyclically shifting the identity matrix.
    
    Parameters:
    Z (int): Dimension of the identity matrix (ZxZ).
    cyclic_value (int): The number of positions to shift; if -1, returns a zero matrix.
    
    Returns:
    numpy.ndarray: The transformed matrix.
    """
    cyclicValue = cyclicValue.astype(int)
    if cyclicValue == -1:
        qout = np.zeros((Z, Z))
    else:
        qout = np.eye(Z)
        qout = np.roll(qout, shift=cyclicValue, axis=1)
    return qout
