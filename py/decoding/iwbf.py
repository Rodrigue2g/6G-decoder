import numpy as np

def ldpc_iwbf(H, re, alp, max_ite):
    # input
    #   H - parity check matrix
    #   re - received word
    #   alp - weighting factor for the weighted bit flipping
    #   max_ite - maximum iteration
    
    # output
    #   syn - syndrome bits
    #   hard - hard-coded message

    # Step 1: Generate the syndrome bits
    row, col = H.shape
    hard = np.zeros(col)
    y_soft = re.copy()    # save soft-decision data for IWBF iteration
    y_re = re.copy()      # save soft-decision data to convert to hard-decision
    iteration = 0

    # hard decision from BPSK
    for i in range(col):
        if y_re[i] > 0.0:
            hard[i] = 1
        else:
            hard[i] = 0
    hard_0 = hard.copy() # save initial hard-decision data

    H_soft = H * np.abs(y_re)

    # Calculate the initial syndrome
    syn = np.mod(hard @ H.T, 2)

    # Step 2: Solve for y_min
    y_min = np.full(row, np.inf) # IWBF y_min per check node/row
    for s1 in range(row):
        for s2 in range(col):
            if H[s1, s2] != 0:
                y_min[s1] = min(y_min[s1], H_soft[s1, s2])

    # Steps 3, 4, and 5: Iterative decoding
    while np.any(syn != 0) and (iteration < max_ite):
        iteration += 1

        # Step 3: Compute for En, IWBF
        En = np.zeros(col)
        for s2 in range(col):
            for s1 in range(row):
                En[s2] += (2*syn[s1]-1)*y_min[s1]*H[s1, s2] - alp*abs(y_re[s2])

        # Step 4: Get maximum En
        id = np.argmax(np.abs(En))

        # Step 5: Flip the bits
        hard[id] = 1 - hard[id]

        # Recompute syndrome bits
        syn = np.mod(hard @ H.T, 2)

    # Check decoding success
    if np.all(syn == 0):
        print('IWBF DECODING IS SUCCESSFUL')
    else:
        print('IWBF DECODING IS UNSUCCESSFUL')

    return syn, hard
