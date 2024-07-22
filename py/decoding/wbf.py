import numpy as np

def sc_ldpc_wbf_decoding(H, chOut, winSize, cpd_L, cpd_w, Z, itera, nb, mb, Nb, Mb, alpha, frames):

    print("H matrix columns (bits):", H.shape[1])
    print("chOut columns (bits per frame):", chOut.shape[0])
    decode_bits = np.zeros((cpd_L * Nb, frames))
    decode_itera = 0

    for i_frame in range(frames):
        hard_decision = np.where(chOut[:, i_frame] > 0, 1, 0)
        if hard_decision.shape[0] != H.shape[1]:
            raise ValueError("Dimension mismatch: hard decision vector length and number of columns in H must match.")
        
        syn = np.mod(np.dot(hard_decision, H.T), 2)
        
        iteration = 0
        while np.any(syn) and iteration < itera:
            iteration += 1
            E = np.zeros_like(chOut[:, i_frame])

            # Iterate over each check node
            for j in range(H.shape[0]):  # for each check
                connected_bits = H[j, :] == 1
                min_value = np.min(np.abs(chOut[connected_bits, i_frame]))
                syn_contribution = (2 * syn[j] - 1) * min_value

                # Update E for all bits connected to check node j
                E[connected_bits] += syn_contribution

            # Determine the bit with the maximum reliability to flip
            flip_bit_index = np.argmax(np.abs(E))
            hard_decision[flip_bit_index] = 1 - hard_decision[flip_bit_index]

            # Update the syndrome
            syn = np.mod(np.dot(hard_decision, H.T), 2)

        decode_bits[:, i_frame] = hard_decision
        decode_itera += iteration

    return decode_bits, decode_itera
