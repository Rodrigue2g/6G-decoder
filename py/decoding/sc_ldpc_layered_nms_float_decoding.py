import numpy as np

def rotation(yIn, shiftNum, Z):
    shift = int(shiftNum)  # Ensure shiftNum is an integer
    return np.roll(yIn, -shift)

def signOP(x):
    return np.sign(x)

def sc_ldpc_layered_nms_float_decoding(chOut, BGWin, BGTerm, winSize, cpd_L, cpd_w, Z, itera, nb, mb, Nb, Mb, alpha, frames):
    BGWin_cols = nb * (winSize + cpd_w - 1)
    BGWin_rows = mb * winSize
    print(f"BGWin_cols {BGWin_cols}")
    print(f"BGWin_rows {BGWin_rows}")
    winDelay = winSize - 1
    run_len = cpd_L + winDelay

    y_full = np.vstack([chOut, np.full((winDelay * Nb - (cpd_w - 1) * Mb, frames), np.inf)])
    decode_bits = np.zeros((cpd_L * Nb, frames))
    decode_itera = 0

    for i_frame in range(frames):
        TermBGCnt = 0
        posteriorMem = np.full((Z, BGWin_cols), 1000)
        extrinsicMem = np.zeros((Z, BGWin_rows, BGWin_cols))

        for i_run in range(run_len):
            posteriorMem[:, :-nb] = posteriorMem[:, nb:]
            extrinsicMem[:, :-mb, :-nb] = extrinsicMem[:, mb:, nb:]

            posteriorMem[:, -nb:] = np.reshape(y_full[(i_run * Nb):((i_run + 1) * Nb), i_frame], (Z, nb))
            extrinsicMem[:, -mb:, :] = 0
            extrinsicMem[:, :, -nb:] = 0

            if i_run > cpd_L:
                BGWD = BGTerm[TermBGCnt * mb: (TermBGCnt + 1) * mb, TermBGCnt * nb: (TermBGCnt + 1) * nb]
                TermBGCnt += 1
            else:
                BGWD = BGWin

            QMem = posteriorMem.copy()
            RMem = extrinsicMem.copy()
            TMem = np.zeros_like(posteriorMem)
            QSign = np.zeros_like(posteriorMem, dtype=int)

            for i_itera in range(itera):
                TDFlag = True
                SynCheck = np.zeros((Z,))

                for i_rows in range(BGWin_rows):
                    min1 = np.full((Z,), np.inf)
                    min2 = np.full((Z,), np.inf)
                    sign = np.ones((Z,))
                    minIdx = np.zeros((Z,), dtype=int)
                    minVec = np.zeros((Z,))

                    for i_cols in range(BGWin_cols):
                        """
                        print(f"i_rows {i_rows}")
                        print(f"i_cols {i_cols}")
                        print(f"BGWD.shape {BGWD.shape}")
                        
                        if i_rows >= BGWD.shape[0]:
                            continue
                        if i_cols >= BGWD.shape[1]:
                            continue
                        """
                        if i_rows < BGWD.shape[0] and i_cols < BGWD.shape[1]:
                            if BGWD[i_rows, i_cols] == -1:
                                continue
                            shiftNum = int(BGWD[i_rows, i_cols])
                            Qc = rotation(QMem[:, i_cols], shiftNum, Z)
                            TMem[:, i_cols] = Qc - RMem[:, i_rows, i_cols]
                            sign *= signOP(TMem[:, i_cols])
                            TMemAbs = np.abs(TMem[:, i_cols])
                            QSign[:, i_cols] = signOP(Qc)

                            for i_Z in range(Z):
                                if TMemAbs[i_Z] < min1[i_Z]:
                                    min2[i_Z] = min1[i_Z]
                                    min1[i_Z] = TMemAbs[i_Z]
                                    minIdx[i_Z] = i_cols
                                elif TMemAbs[i_Z] < min2[i_Z]:
                                    min2[i_Z] = TMemAbs[i_Z]

                    for i_cols in range(BGWin_cols):
                        if i_rows < BGWD.shape[0] and i_cols < BGWD.shape[1]:
                            if BGWD[i_rows, i_cols] == -1:
                                continue
                            shiftNum = int(BGWD[i_rows, i_cols])
                            if i_cols < (cpd_w - 1) * nb:
                                Qtmp = rotation(QMem[:, i_cols], shiftNum, Z)
                                SynCheck = np.logical_xor(SynCheck, (Qtmp < 0))
                            else:
                                for i_Z in range(Z):
                                    minVec[i_Z] = min2[i_Z] if minIdx[i_Z] == i_cols else min1[i_Z]
                                RMem[:, i_rows, i_cols] = sign * signOP(TMem[:, i_cols]) * minVec * alpha
                                Qtmp = TMem[:, i_cols] + RMem[:, i_rows, i_cols]
                                if np.any((Qtmp < 0) != (QSign[:, i_cols] < 0)):
                                    TDFlag = False
                                SynCheck = np.logical_xor(SynCheck, (Qtmp < 0))
                                QMem[:, i_cols] = rotation(Qtmp, Z - shiftNum, Z)

                    if np.sum(SynCheck) > 0:
                        TDFlag = False

                if TDFlag:
                    break

            posteriorMem = QMem
            extrinsicMem = RMem

            if i_run > winDelay:
                windowOut = (posteriorMem < 0).astype(float).reshape(BGWin_cols * Z)
                decode_bits[(i_run - winDelay - 1) * Nb: (i_run - winDelay) * Nb, i_frame] = windowOut[(cpd_w - 1) * Nb: cpd_w * Nb]
                decode_itera += i_itera

    return decode_bits, decode_itera
