
import numpy as np

def gdbf_decoder(H_full, input_data, L, w, Z, nb, mb, itera=60, tail_length=128):
    #SC-LDPC full matrix dimensions
    M_tot = mb * Z * L + tail_length
    N_tot = nb * Z * L + tail_length
    #LDPC H matrix dimensions
    M = mb * Z
    N = w * nb * Z

    code_length = nb * Z

    H = np.zeros((M, N))
    H_full_flat = np.copy(H_full)#.flatten()
    #print(f"H_full_flat.shape {H_full_flat.shape}\n")
    """
    for i in range(M):
        for j in range(N):
            H[i][j] = H_full_flat[(M+i)*N_tot+j]
            #H[i][j] = H_full[(M+i)][N_tot+j]
    """
    H = H_full_flat[M : 2*M , 0 : N]
    print(f"H {H[10][1]}")
    

    # Convert input to +1/-1
    xhat_tot = np.ones(N_tot+code_length, dtype=int)
    for i in range(N_tot):
        #xhat_tot[i] = (input_data[i]>=0 ? 1 : -1)
        if input_data[i]>=0:
            xhat_tot[i] = 1
        elif input_data[i]<0:
            xhat_tot[i] = -1

    
    # Initialize the output array
    output = np.zeros(input_data.shape)
    #output = np.ones(N_tot+code_length, dtype=int)
    # Process each codeword
    for t in range(L-1):
        """Spatialy coupled"""
        y = input_data[t*code_length :]
        xhat = xhat_tot[t*code_length :]
        
        #print(f"xhat.shape {xhat.shape}")
        for i_itera in range(itera):
            # Check parity
            holds = True
            for i in range(M):
                accum = 1
                for j in range(N):
                    if(H[i][j] == 1):
                        accum *= xhat[j]
                if(accum != 1):
                    holds = False
                    break
            if holds:
                print(f"Converged after {i_itera} iterations")
                break
            elif i_itera == itera - 1:
                print("Did not converge after maximum iterations")

            # Gradient descent to adjust xhat
            lowest_value = 0 #float('inf')
            lowest_idx = 0
            for k in range(N):
                inver_func = xhat[k] * y[k]
                for mk in range(M):
                    if H[mk][k] == 1:
                        accum = 1
                        for nk in range(N):
                            if H[mk][nk] == 1:
                                accum *= xhat[nk]
                        inver_func += accum

                # Find the smallest value in inverFunc
                if (inver_func < lowest_value) or (k==0):
                    lowest_idx = k
                    lowest_value = inver_func

            # Flip the bit with the smallest value
            xhat[lowest_idx] = -xhat[lowest_idx]

        # Update output decoded bits
        for i in range(code_length):
            output[t*code_length+i] = (xhat[i] < 0)

    return output, itera

