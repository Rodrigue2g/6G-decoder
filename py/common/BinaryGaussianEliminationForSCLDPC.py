import numpy as np

def BinaryGaussianEliminationForSCLDPC(H, M, N, w):
    K = N - M
    Identity = np.eye(M)
    IdentityInv = np.eye(M)
    
    # Move last Mb columns to the beginning
    H[:, :, 0] = np.hstack([H[:, K:, 0], H[:, :K, 0]])
    
    if np.array_equal(H[:, :M, 0], np.eye(M)):
        print('Original matrix has already satisfied a standard format.')
    else:
        print('Perform Gaussian elimination by row transformation to get a standard format.')

        # Start to do Gaussian Elimination
        for i in range(M):
            pivot = -1
            for j in range(i, M):
                if H[j, i, 0] == 1:
                    pivot = j
                    break
            
            if pivot != -1:
                # Add
                if i != pivot:
                    H[i, :, :w] = (H[i, :, :w] + H[pivot, :, :w]) % 2
                    Identity[i, :] = (Identity[i, :] + Identity[pivot, :]) % 2
                    IdentityInv[:, pivot] = (IdentityInv[:, pivot] + IdentityInv[:, i]) % 2
                
                # Eliminate other rows (from top to bottom)
                for j in range(i + 1, M):
                    if H[j, i, 0] == 1:
                        H[j, :, :w] = (H[j, :, :w] + H[i, :, :w]) % 2
                        Identity[j, :] = (Identity[j, :] + Identity[i, :]) % 2
                        IdentityInv[:, i] = (IdentityInv[:, i] + IdentityInv[:, j]) % 2
            else:
                print(f"Cannot find a valid pivot for row {i}")

        # Eliminate other rows (from bottom to top)
        for i in range(M-1, -1, -1):
            for j in range(i-1, -1, -1):
                if H[j, i, 0] == 1:
                    H[j, :, :w] = (H[j, :, :w] + H[i, :, :w]) % 2
                    Identity[j, :] = (Identity[j, :] + Identity[i, :]) % 2
                    IdentityInv[:, i] = (IdentityInv[:, i] + IdentityInv[:, j]) % 2

    # Move back
    H[:, :, 0] = np.hstack([H[:, M:, 0], H[:, :M, 0]])
    return H, Identity, IdentityInv
