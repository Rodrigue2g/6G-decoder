import numpy as np
from scipy.io import loadmat
from numpy.linalg import inv
import matplotlib.pyplot as plt
import time

from sc_ldpc_layered_nms_float_decoding import sc_ldpc_layered_nms_float_decoding

import sys
sys.path.append('common')
from BaseGraphTermAuxiliaryMatrixGenerator import BaseGraphTermAuxiliaryMatrixGenerator
from BinaryGaussianEliminationForSCLDPC import BinaryGaussianEliminationForSCLDPC
from identityMatrixExtend import identityMatrixExtend
from signOP import signOP

start_time = time.time()
np.random.seed(1)

# Parameter definition
cpd_L = 50
cpd_w = 2  # also defined as ms or mcc = w-1
winSize = 5
itera = 10

nb = 40
mb = 8
Z = 16
kb = nb - mb

Nb = nb * Z
Mb = mb * Z
Kb = kb * Z

frames = 500
groups = 10
MultiCoreFlag = 1
max_frames = int(1e+9)
max_errors = 80

alpha = 0.75
beta = 0.5

snr = np.arange(2, 3.375, 0.125)
rate = cpd_L * Kb / (cpd_L * Nb + (cpd_w-1) * Mb)

nb_efer = np.zeros(len(snr))
nb_ebler = np.zeros(len(snr))
nb_eber = np.zeros(len(snr))
nb_itera = np.zeros(len(snr))
nb_frame = np.zeros(len(snr))

FER = np.zeros(len(snr))
BLER = np.zeros(len(snr))
BER = np.zeros(len(snr))
ITER = np.zeros(len(snr))

framesNum = 0

LoadFullPath = f'construction/RPTU_SCLDPC_Mb{mb}_Nb{nb}_Z{Z}.mat'
profile = loadmat(LoadFullPath)
BGCpd = np.zeros((mb, nb, cpd_w))  # construct basic matrices
BGCpd[:, :, 0] = profile['BGCpd'][:, :, 0]
BGCpd[:, :, 1] = profile['BGCpd'][:, :, 1]

#generate the coupled base graph for window
BGWin_len = winSize + cpd_w - 1
BGWin = np.zeros((mb * winSize, nb * BGWin_len)) - 1
#BGWin = np.full((mb * winSize, nb * BGWin_len), -1)
"""
print(f"BGWin_len: {BGWin_len}")
print(f"BGWin shape: {BGWin.shape}")
print(f"BGCpd shape: {BGCpd.shape}")
"""
for i_row in range(1, winSize+1):
    for i_w in range(1, cpd_w+1):
        #BGWin((i_row-1)*mb+1 : i_row*mb, (i_w+i_row-1-1)*nb+1 : (i_w+i_row-1)*nb) = BGCpd(:, :, cpd_w-i_w+1);
        """
        print(f" ")
        print(f"1 - (i_row-1)*mb+1: {(i_row-1)*mb+1}")
        print(f"2 - i_row*mb: {i_row*mb}")

        print(f"3 - (i_w+i_row-1-1)*nb+1: {(i_w+i_row-1-1)*nb+1}")
        print(f"4 -  (i_w+i_row-1)*nb: {(i_w+i_row-1)*nb}")

        print(f"5 - cpd_w - i_w+1: {cpd_w - i_w+1}")
        """
        if i_row != winSize:
            BGWin[(i_row-1)*mb+1 : i_row*mb + 1, (i_w+i_row-1-1)*nb  : (i_w+i_row-1)*nb] = BGCpd[:, :, cpd_w - i_w]
        elif i_row == winSize:
            """
            Weird as this should be from 33 to 40 but it only gives the right shape from 32 to 40
            print(f"1 - (i_row-1)*mb+1: {(i_row-1)*mb}")
            print(f"2 - i_row*mb: {i_row*mb}")
            print(f"BGWin shape: {BGWin[(i_row-1)*mb+1 : i_row*mb, (i_w+i_row-1-1)*nb  : (i_w+i_row-1)*nb].shape}")
            """
            BGWin[(i_row-1)*mb : i_row*mb, (i_w+i_row-1-1)*nb  : (i_w+i_row-1)*nb] = BGCpd[:, :, cpd_w - i_w]

"""
print(f"BGWin_len: {BGWin_len}")
print(f"BGWin shape: {BGWin.shape}")
print(f"BGCpd shape: {BGCpd.shape}")
print("BGWin content:")
print(BGWin[:,:])
"""

# Generate an auxiliary BG for termination
BGTerm = BaseGraphTermAuxiliaryMatrixGenerator(BGCpd, winSize, cpd_w, BGWin_len, mb, nb)

HCpd = np.zeros((Mb, Nb, cpd_w))
for i_w in range(cpd_w):
    for ii in range(mb):
        for jj in range(nb):
            HCpd[(ii*Z):(ii+1)*Z, (jj*Z):(jj+1)*Z, i_w] = identityMatrixExtend(Z, BGCpd[ii, jj, i_w])

# Modify H0 as the standard format
HCpd, IdentityTail, IdentityTailInv = BinaryGaussianEliminationForSCLDPC(HCpd, Mb, Nb, cpd_w)

# Extract H0'
H0Prime = HCpd[:, :Kb, 0]


for i_run in range(0, max_frames, frames):
    if nb_efer[-1] >= max_errors and nb_ebler[-1] >= 100:
        break

    # Generate random binary input data
    u = np.random.randint(0, 2, (Kb, frames, cpd_L))
    x = np.zeros((Nb, frames, cpd_L))
    partialSyndrome = np.zeros((Mb, frames, cpd_w-1))

    # Encoding process
    for i_L in range(cpd_L):
        p = (np.dot(H0Prime, u[:, :, i_L]) + partialSyndrome[:, :, 0]) % 2
        x[:, :, i_L] = np.vstack([u[:, :, i_L], p])

        # Implement cyclic shifting registers for partial syndromes
        partialSyndrome_last = np.copy(partialSyndrome)
        partialSyndrome_block = np.zeros((Mb, frames, cpd_w-1))
        
        for i_w in range(cpd_w - 1):
            partialSyndrome_block[:, :, i_w] = (np.dot(HCpd[:, :, i_w + 1], x[:, :, i_L])) % 2

        # Cyclical shifting and updating
        for i_w in range(cpd_w - 2):
            partialSyndrome[:, :, i_w] = (partialSyndrome_last[:, :, i_w + 1] + partialSyndrome_block[:, :, i_w]) % 2
        partialSyndrome[:, :, -1] = partialSyndrome_block[:, :, -1]

    # Permute and reshape for transmission
    x_sent_full = np.transpose(x, (0, 2, 1)).reshape(cpd_L * Nb, frames)
    
    # Add a short tail composed of partial syndrome
    for i_w in range(cpd_w-1):
        partialSyndrome[:, :, i_w] = (np.dot(IdentityTailInv, partialSyndrome[:, :, i_w])) % 2
        x_sent_full = np.vstack([x_sent_full, partialSyndrome[:, :, i_w]])

    # Noise generation (AWGN)
    noise = np.random.normal(0, 1, x_sent_full.shape)
    for i_snr, db in enumerate(snr):
        if nb_efer[i_snr] >= max_errors and nb_ebler[i_snr] >= 100:
            continue

        # Channel: AWGN and BPSK modulation
        sigma = 1 / np.sqrt(2.0 * rate) * 10 ** (-db / 20)
        symbol = 1 - 2 * x_sent_full
        waveform = symbol + noise * sigma
        chOut = 2 * waveform / sigma**2

        # Avoid the risk of array out of bounds in input LLR vector
        chOut_full = np.vstack([chOut, 1000 * np.ones(((winSize - 1) * Nb - (cpd_w-1) * Mb, frames))])

        # Decoding function
        # [rxcbs, rxitera] = sc_ldpc_layered_nms_float_decoding_mex(chOut_full, BGWin, BGTerm, winSize, cpd_L, cpd_w, Z, itera, nb, mb, Nb, Mb, alpha, frames);
        
        rxcbs, rxitera = sc_ldpc_layered_nms_float_decoding(chOut_full, BGWin, BGTerm, winSize, cpd_L, cpd_w, Z, itera, nb, mb, Nb, Mb, alpha, frames)

        # Logic to calculate BER, FER, ...
        # Initialization of matrices to store decoded bits and original bits for comparison
        dec_info_bits = np.zeros((Kb, cpd_L, frames))
        x_info_bits = np.zeros((Kb, cpd_L, frames))

        # Extracting decoded and original bits for comparison
        for i_L in range(cpd_L - 1):
            dec_info_bits[:, i_L, :] = rxcbs[(i_L * Nb):(i_L * Nb + Kb), :]
            x_info_bits[:, i_L, :] = x_sent_full[(i_L * Nb):(i_L * Nb + Kb), :]

        # Calculate number of erroneous bits
        errNum = np.sum(np.abs(x_info_bits - dec_info_bits), axis=2)

        # Update error metrics
        nb_eber[i_snr] += np.sum(errNum)
        nb_ebler[i_snr] += np.sum(errNum > 0)
        nb_efer[i_snr] += np.sum(np.any(errNum > 0, axis=0))
        nb_itera[i_snr] += rxitera
        nb_frame[i_snr] += frames

        # Calculate rates
        valid_indices = nb_frame > 0  # array of boolean

        # Fill with NaN for invalid computations
        FER = np.full_like(nb_frame, np.nan, dtype=float)
        BLER = np.full_like(nb_frame, np.nan, dtype=float)
        BER = np.full_like(nb_frame, np.nan, dtype=float)
        ITER = np.full_like(nb_frame, np.nan, dtype=float)

        FER[valid_indices] = nb_efer[valid_indices] / nb_frame[valid_indices]
        BLER[valid_indices] = nb_ebler[valid_indices] / (nb_frame[valid_indices] * cpd_L)
        BER[valid_indices] = nb_eber[valid_indices] / (nb_frame[valid_indices] * cpd_L * Kb)
        ITER[valid_indices] = nb_itera[valid_indices] / (nb_frame[valid_indices] * cpd_L)

        framesNum += frames

        # Performance information
        print(f"\nRunning {framesNum} frames")
        print("Current block error performance:")
        print(f"SNR: {snr}, FER: {FER}, BLER: {BLER}, BER: {BER}, ITER: {ITER}")
        #print(f"SNR: {db:.2f} dB, FER: {FER:.5f}, BLER: {BLER:.5f}, BER: {BER:.5f}, ITER: {ITER:.2f}")

plt.figure(figsize=(13, 6))
plt.semilogy(snr, FER, 'o-b', linewidth=1)
plt.title('FER')
plt.grid(True)

plt.figure(figsize=(13, 6))
plt.semilogy(snr, BLER, 'o-b', linewidth=1)
plt.title('BLER')
plt.grid(True)

plt.figure(figsize=(13, 6))
plt.semilogy(snr, BER, 'o-r', linewidth=1)
plt.xlabel('SNR', fontsize=14)
plt.ylabel('BER', fontsize=14)
plt.title('BER vs SNR')
plt.grid(True, linestyle='--')
plt.show()

print("Total runtime:", time.time() - start_time)
