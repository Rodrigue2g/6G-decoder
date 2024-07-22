import numpy as np

def hybrid_weighted_bit_flipping(H, r, max_iter=10):
    # Convert received vector r into hard decisions (0 or 1)
    r_hard = (r > 0).astype(int)  # Assuming BPSK: 1 -> 0, -1 -> 1
    
    # Initialize the codeword estimate with the hard decisions
    c = r_hard.copy()
    
    # Calculate the number of checks each bit is involved in (weight for flipping)
    weights = np.sum(H, axis=0)
    
    # Main decoding loop
    for _ in range(max_iter):
        # Compute the syndrome
        s = H.dot(c) % 2
        print(f"syndrome {s}\n")

        # Check if all parities are satisfied
        if np.all(s == 0):
            break
        
        # Compute the error metrics for each variable node
        e = np.zeros_like(c)
        for j in range(H.shape[1]):  # For each variable node
            for i in range(H.shape[0]):  # For each check node
                if H[i, j] == 1:
                    # Contribution of each check node
                    e[j] += s[i] * weights[j]
        
        # Flip the bits with the highest error metrics
        flip_threshold = np.percentile(e, 90)  # Flipping threshold, e.g., top 10%
        c[e >= flip_threshold] = 1 - c[e >= flip_threshold]
    
    return c