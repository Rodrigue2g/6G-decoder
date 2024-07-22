import numpy as np

"""
Decoding Steps of the Proposed Low Complex 
Hybrid Weighted Bit Flipping Algorithm
"""
def lc_hwbf_decoder(H_full, input_data, L, w, Z, nb, mb, max_iter=10, tail_length=128):
    """
    Implements a Hybrid Weighted Bit-Flipping algorithm for SC-LDPC codes.

    Parameters:
        H (numpy array): Parity check matrix.
        y (numpy array): Received LLRs.
        L (int): Number of spatial positions in the SC-LDPC code.
        w (int): Coupling width.
        Z (int): Expansion factor.
        nb (int): Number of columns in the base matrix.
        mb (int): Number of rows in the base matrix.
        max_iter (int): Maximum number of decoding iterations.
        tail_length (int): Length of the tail in the code.

    Returns:
        numpy array: Decoded binary vector.
        int: Number of iterations executed.
    """
    M_tot = mb * Z * L
    N_tot = nb * Z * L
    #LDPC H matrix dimensions
    M = mb * Z
    N = w * nb * Z

    code_length = nb * Z

    """Initialisation of H & xhat"""
    H = np.zeros((M, N))
    H_full_flat = H_full.flatten(order='F')
    for i in range(M):
        for j in range(N):
            H[i][j] = H_full_flat[(M+i)*N_tot+j]

    # Convert input to +1/-1
    xhat_tot = np.ones(N_tot+code_length, dtype=int)
    for i in range(0, N_tot):
        if input_data[i]>=0:
            xhat_tot[i] = 1
        elif input_data[i]<0:
            xhat_tot[i] = -1

    """
    for i in range(M-1):
        if i not in n_indices:
            w_mn[i] = np.min(np.abs(y[i]))
    """

    """Step 2"""
    k = 0
    k_max = 20
    output = np.zeros(input_data.shape)
    Zk = np.zeros(input_data.shape)
    for k in range(k_max):
        """Shift accordingly to the spatial position"""
        y = input_data[k*code_length :]
        out = output[k*code_length :]
        xhat = xhat_tot[k*code_length :]
        # Initial hard decision
        Zk[k*code_length :] = (y < 0).astype(int)

        """Step 1"""
        # Precompute min values for each check node
        w_m = np.zeros(M)
        w_mn = np.zeros((M, N))
        for m in range(M-1):
            # n_indices : indices of min{|y_n|}
            n_indices = np.where(H[m, :] == 1)[0]
            w_m[m] = np.min(np.abs(y[n_indices]))
            for n in n_indices:
                w_mn[m, n] = np.min(np.abs(y[np.setdiff1d(n_indices, n)]))
            
        """Step 3"""
        Sm = Zk @ H % 2
        if np.all(Sm == 0):
            break 
        """Step 4"""
        e = np.zeros(N_tot)
        for n in range(N_tot-1):
            connected_checks = np.where(H[:, n] == 1)[0]
            e[n] = sum((2 * Sm[connected_checks] - 1) * (np.min(np.abs(y[connected_checks])))) / np.abs(y[n])
            """
            m_indices = np.where(H[:, n] == 1)[0]
            e[n] = 1 / np.abs(y[n]) * sum((2 * Sm[connected_checks] - 1) * w_mn[connected_checks, n] * theta_attn)
            """

        """Step 5"""
        Zk = np.where(e <= 0, 0, 1)

        """Step 6"""
        n_flip = np.argmax(np.abs(e))
        Zk[n_flip] = 1 - Z[n_flip]
    
    return Zk, k