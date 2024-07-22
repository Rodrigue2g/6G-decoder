% ---------------------------------------------------------------------------------
% FUNCTION INFORMATION (c) 2023 Telecommunications Circuits Laboratory, EPFL
% ---------------------------------------------------------------------------------
% name  : sc_ldpc_layered_nms_float_decoding
% descr : decoding of spatially coupled LDPC codes

function [decode_bits, decode_itera] = sc_ldpc_layered_nms_float_decoding(chOut, BGWin, BGTerm, winSize, cpd_L, cpd_w, Z, itera, nb, mb, Nb, Mb, alpha, frames)
% ---------------------------------------------------------------------------------
% parameter definition of window
% ---------------------------------------------------------------------------------
BGWin_cols = nb*(winSize + cpd_w - 1);
BGWin_rows = mb*winSize;

% note that this delay is much important in SC-LDPC decoding, which means the latency of sliding in and sliding out
% part1 : due to the latency of loading window, we can not output at first delay steps
% part2 : we need special processing at last delay steps to deal with terminated bits
winDelay = winSize - 1;
run_len  = cpd_L + winDelay; % running steps

% complement the input LLR vector to avoid the risk of array out of bounds
y_full = [chOut; ones(winDelay*Nb - (cpd_w-1)*Mb, frames)*Inf];

% decoded output
decode_bits  = zeros(cpd_L*Nb, frames);
decode_itera = 0;

for i_frame = 1 : frames
    TermBGCnt = 0;
    
    % posterior and extrinsic memory
    posteriorMem = ones(Z, BGWin_cols)*1000; % Note that we CANNOT set as Inf, which may conflict with the following minima
    extrinsicMem = zeros(Z, BGWin_rows, BGWin_cols);

    for i_run = 1 : run_len
        % shift two memories
        posteriorMem(:, 1 : end-nb)             = posteriorMem(:, nb+1 : end);             % from right to left
        extrinsicMem(:, 1 : end-mb, 1 : end-nb) = extrinsicMem(:, mb+1 : end, nb+1 : end); % from down right to upper left

        % load channel LLRs
        posteriorMem(:, end-nb+1 : end) = reshape(y_full((i_run-1)*Nb+1 : i_run*Nb, i_frame), Z, nb); % load channel LLRs
        %extrinsicMem(:, end-mb+1 : end, end-nb+1 : end) = 0; % clean
        extrinsicMem(:, end-mb+1 : end, :) = 0; % clean
        extrinsicMem(:, :, end-nb+1 : end) = 0; % clean

        % load PCM of sliding window
        if i_run > cpd_L
            BGWD      = BGTerm(1 + TermBGCnt*mb : BGWin_rows + TermBGCnt*mb, 1 + TermBGCnt*nb : BGWin_cols + TermBGCnt*nb);
            TermBGCnt = TermBGCnt + 1;
        else
            BGWD = BGWin;
        end

        % start decoding
        QMem = posteriorMem;
        RMem = extrinsicMem;

        % temporary memory
        TMem  = zeros(Z, BGWin_cols);
        QSign = zeros(Z, BGWin_cols);

        for i_itera = 1 : itera
            TDFlag = 1; % partial syndrome check

            for i_rows = 1 : BGWin_rows

                SynCheck = zeros(Z, 1);

                % -------------------------------------------------------------------
                % Phase 1: MIN
                % -------------------------------------------------------------------
                min1   = Inf * ones(Z, 1); % store the first minima
                min2   = Inf * ones(Z, 1); % store the second minima
                sign   = ones(Z, 1);       % store the sign bit of TMem (Z-based)

                minIdx = zeros(Z, 1);
                minVec = zeros(Z, 1);

                for i_cols = 1 : BGWin_cols
                    % checkout whether there is -1 in base matrix, if yes, directly skip
                    if BGWD(i_rows, i_cols) == -1
                        continue;
                    end
                    shiftNum = BGWD(i_rows, i_cols); % cyclic shifting value (c of P^c)

                    % Calculate the Qc and TMem
                    Qc = Rotation(QMem(:, i_cols), shiftNum, Z);    % rotate left by shiftNum bits
                    TMem(:, i_cols) = Qc - RMem(:, i_rows, i_cols); % TMem is actually the previous Qvc (from VN --> CN)

                    % Get the sign and magnitude
                    sign             = sign.*signOP(TMem(:, i_cols));
                    TMemAbs          = abs(TMem(:, i_cols));
                    QSign(:, i_cols) = signOP(Qc);

                    % Search the min1, min2, and the column index corresponding to min1
                    for i_Z = 1 : Z
                        if TMemAbs(i_Z) < min1(i_Z)
                            min2(i_Z)   = min1(i_Z);
                            min1(i_Z)   = TMemAbs(i_Z);
                            minIdx(i_Z) = i_cols;
                        elseif TMemAbs(i_Z) < min2(i_Z)
                            min2(i_Z)   = TMemAbs(i_Z);
                        end
                    end
                end

                % -------------------------------------------------------------------
                % Phase 2: Q- and R- messages update
                % -------------------------------------------------------------------
                for i_cols = 1 : BGWin_cols
                    % checkout whether there is -1 in base matrix, if yes, directly skip
                    if BGWD(i_rows, i_cols) == -1
                        continue;
                    end
                    shiftNum = BGWD(i_rows, i_cols); % cyclic shifting value (c of P^c)

                    % work on the previous ms decoded blocks, just execute parity check w.o. updating messages
                    if i_cols <= (cpd_w-1)*nb
                        Qtmp = Rotation(QMem(:, i_cols), shiftNum, Z);

                        SynCheck = xor(SynCheck, (Qtmp < 0));
                    else
                        % Deal with the special cases that the minina belongs to the current column
                        for i_Z = 1 : Z
                            if minIdx(i_Z) == i_cols
                                minVec(i_Z) = min2(i_Z);
                            else
                                minVec(i_Z) = min1(i_Z);
                            end
                        end

                        RMem(:, i_rows, i_cols) = sign.*signOP(TMem(:, i_cols)).*minVec*alpha;
                        Qtmp = TMem(:, i_cols) + RMem(:, i_rows, i_cols);

                        if sum(abs((Qtmp < 0) - (QSign(:, i_cols) < 0))) ~= 0
                            TDFlag = 0;
                        end
                        SynCheck = xor(SynCheck, (Qtmp <0));

                        % Inverse rotation, i.e., P^c*P^c' = I
                        Qtmp = Rotation(Qtmp, Z - shiftNum, Z);
                        QMem(:, i_cols) = Qtmp;
                    end
                end
                
                if sum(SynCheck) > 0
                    TDFlag = 0;
                end

            end

            if TDFlag == 1
                break;
            end
        end

        posteriorMem = QMem;
        extrinsicMem = RMem;

        if i_run > winDelay
            windowOut = double(reshape(posteriorMem, BGWin_cols*Z, 1) < 0);
            decode_bits((i_run - winDelay - 1)*Nb + 1 : (i_run - winDelay)*Nb, i_frame) = windowOut((cpd_w - 1)*Nb + 1 : cpd_w*Nb);
            decode_itera = decode_itera + i_itera;
        end
    end
end
end

function yOut = Rotation(yIn, shiftNum, Z)
yOut = [yIn(shiftNum + 1 : Z); yIn(1 : shiftNum)];
end