import numpy as np

def hard_decoder(H, c, MAX_ITER):
    """
    Parameters:
        c (array): Binary column vector [1, N] -> input codeword
        H (array): Binary matrix [M, N] (True and False)
        MAX_ITER (int): Maximum number of iterations
    
    Returns:
        c_cor (array): Binary column vector [1, N] resulting from decoding
    """
    
    # Transpose inputs since Python uses different array orientation
    c = np.array(c).reshape(-1, 1)  # Ensure c is a column vector
    
    nCheckNodes = H.shape[0]  # Number of check nodes (rows in H)
    nVariableNodes = H.shape[1]  # Number of variable nodes (columns in H)
    
    # Parity matrix for parity checks
    parity = np.zeros(nCheckNodes)
    c_cor = c.flatten()  # Output vector
    nIter = 1  # Iteration count
    
    # Calculate the majority threshold for each variable node
    majority = np.zeros(nVariableNodes)
    for j in range(nVariableNodes):
        majority[j] = (np.sum(H[:, j]) + 1) // 2
    
    while nIter <= MAX_ITER: # and np.sum(c_cor) % 2 == 1:
        print(f"c_cor[j] {c_cor[j]}\n")
        # While maximum iterations not exceeded and parity test fails
        for i in range(nCheckNodes):
            parity[i] = np.mod(np.dot(H[i, :], c_cor), 2)  # Parity test for each check node
        
        for j in range(nVariableNodes):
            nOnes = c_cor[j]  # Initial count of ones for the variable node
            for i in range(nCheckNodes):
                if H[i, j] == 1:  # If the check node is connected to the variable node
                    nOnes += np.mod(parity[i] + c_cor[j], 2)  # Account for the check node's value
            
            # Decide the new value of the variable node based on majority
            if nOnes > majority[j]:
                c_cor[j] = 1
            else:
                c_cor[j] = 0
        
        nIter += 1
    
    # Ensure output is a column vector as required
    #c_cor = c_cor.reshape(-1, 1)
    #for i in range(len(c_cor)):
    #    c_cor[i] = (c_cor[i]<0)
    return c_cor
             
