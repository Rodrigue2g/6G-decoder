import numpy as np

def BaseGraphTermAuxiliaryMatrixGenerator(BaseGraph, winSize, cpd_w, BaseGraphWin_len, mb, nb):
    """
    Generate the auxiliary base graph for the last terminated bits.
    
    Parameters:
    BaseGraph (numpy.ndarray): The base graph matrix of size (mb, nb, cpd_w).
    winSize (int): Window size for processing.
    cpd_w (int): Number of components or windows in the base graph processing.
    BaseGraphWin_len (int): Length of the base graph window.
    mb (int): Number of rows in each component of the base graph.
    nb (int): Number of columns in each component of the base graph.
    
    Returns:
    numpy.ndarray: The auxiliary base graph termination matrix.
    """
    winDelay = winSize - 1  # delay of sliding in and sliding out

    BaseGraphTerm_cols = BaseGraphWin_len + winDelay - 1
    BaseGraphTerm_rows = winSize + winDelay - 1
    BaseGraphTermMatrix = np.full((BaseGraphTerm_rows * mb, BaseGraphTerm_cols * nb), -1)

    # Fill the component matrices
    for i_rows in range(1, BaseGraphTerm_rows+1):
        for i_w in range(1, cpd_w + 1):
            #BaseGraphTermMatrix((i_rows-1)*mb+1 : i_rows*mb, (i_rows+i_w-2)*nb+1 : (i_rows+i_w-1)*nb) = BaseGraph(:, :, cpd_w-i_w+1);
            if i_rows != BaseGraphTerm_rows:
                BaseGraphTermMatrix[(i_rows-1)*mb+1 : i_rows*mb+1, (i_rows+i_w-2)*nb : (i_rows+i_w-1)*nb] = BaseGraph[:, :, cpd_w - i_w]
            elif i_rows == BaseGraphTerm_rows:
                """
                Same problem as in BGWin generation, should go from 57 to 64 but only works from 56 to 65
                otw shape is (7, 40) instead of (8, 40)
                print(f"1 - (i_rows-1)*mb+1: {(i_rows-1)*mb+1}")
                print(f"2 - i_rows*mb: {i_rows*mb}")
                print(f"BaseGraphTermMatrix shape: {BaseGraphTermMatrix[(i_rows-1)*mb : i_rows*mb+1, (i_rows+i_w-2)*nb : (i_rows+i_w-1)*nb].shape}")
                """
                BaseGraphTermMatrix[(i_rows-1)*mb : i_rows*mb+1, (i_rows+i_w-2)*nb : (i_rows+i_w-1)*nb] = BaseGraph[:, :, cpd_w - i_w]

    # Fill the last columns with -1
    BaseGraphTermMatrix[:, (BaseGraphWin_len - 1) * nb:] = -1

    # Add an identity matrix at the designated area
    IdentityTail = np.eye(mb * (cpd_w - 1)) - 1
    IdentityTailStartRow = (winSize - 1) * mb
    IdentityTailEndRow = IdentityTailStartRow + mb * (cpd_w - 1)
    IdentityTailStartCol = (BaseGraphWin_len - 1) * nb
    IdentityTailEndCol = IdentityTailStartCol + (cpd_w - 1) * mb
    BaseGraphTermMatrix[IdentityTailStartRow:IdentityTailEndRow, IdentityTailStartCol:IdentityTailEndCol] = IdentityTail

    return BaseGraphTermMatrix
