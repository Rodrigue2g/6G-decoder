% ---------------------------------------------------------------------------------
% FUNCTION INFORMATION (c) 2023 Telecommunications Circuits Laboratory, EPFL
% ---------------------------------------------------------------------------------
% name  : CouplingSpreadingRandom
% descr : use random coupling spreading to transform 5G-NR LDPC to SC-LDPC

rng(10);
TxRx.bgn   = 1;
TxRx.cpd_w = 2;
TxRx.Z     = 32;
TxRx.Kb    = 22;
TxRx.Mb    = 22;
TxRx.Nb    = TxRx.Mb + TxRx.Kb;

LoadfullPath = sprintf('construction\\ParityCheckMatrixProfiles\\5gNR_BG%d_Mb%d_Nb%d_Z%d.mat', TxRx.bgn, TxRx.Mb, TxRx.Nb, TxRx.Z);
profile  = load(LoadfullPath);

BGCpd = zeros(TxRx.Mb, TxRx.Nb, TxRx.cpd_w) - 1;
randPattern = rand(TxRx.Mb, TxRx.Kb) > 0.5;

for i_rows = 1:TxRx.Mb
    for i_cols = 1:TxRx.Nb
        if i_cols <= TxRx.Kb
            if randPattern(i_rows, i_cols) == 1
                BGCpd(i_rows, i_cols, 1) = profile.BG_2D(i_rows, i_cols);
            else
                BGCpd(i_rows, i_cols, 2) = profile.BG_2D(i_rows, i_cols);
            end
        else
            BGCpd(i_rows, i_cols, 1) = profile.BG_2D(i_rows, i_cols);
        end
    end
end

SaveFullPath  = sprintf('construction\\ParityCheckMatrixProfiles\\RandSpread_5gNR_BG%d_Mb%d_Nb%d_Z%d.mat', TxRx.bgn, TxRx.Mb, TxRx.Nb, TxRx.Z);
save(SaveFullPath, 'BGCpd');