% ---------------------------------------------------------------------------------
% FUNCTION INFORMATION (c) 2023 Telecommunications Circuits Laboratory, EPFL
% ---------------------------------------------------------------------------------
% name  : spatially_coupled_ldpc_platform
% descr : encoding and decoding of time-invarient spatially-coupled LDPC codes

% ---------------------------------------------------------------------------------
% parameter definition
% ---------------------------------------------------------------------------------
tic
rng(1);
cpd_L   = 50;
cpd_w   = 2; % also defined as ms or mcc = w-1
winSize = 5;
itera   = 10;

nb  = 40;
mb  = 8;
Z   = 16;
kb  = nb-mb;

Nb = nb*Z;
Mb = mb*Z;
Kb = kb*Z;

frames        = 1;
groups        = 10;
MultiCoreFlag = 1;
max_frames    = 1000;
max_errors    = 80;

alpha  = 0.75;
beta   = 0.5;

snr   = 2:0.25:8.5; %3.25;

% (cpd_w-1)*Mb is the length of tail we attach an identity matrix 
% at the end to do the termination
rate  = cpd_L*Kb / (cpd_L*Nb + (cpd_w-1)*Mb);

nb_efer  = zeros(1, length(snr));
nb_ebler = zeros(1, length(snr));
nb_eber  = zeros(1, length(snr));
nb_itera = zeros(1, length(snr));
nb_frame = zeros(1, length(snr));

nb_sent_ber =  zeros(1, length(snr));

FER  = zeros(1, length(snr));
BLER = zeros(1, length(snr));
BER  = zeros(1, length(snr));
ITER = zeros(1, length(snr)); 

SENT = zeros(1, length(snr)); 

framesNum = 0;

% ---------------------------------------------------------------------------------
% Hw generation
% ---------------------------------------------------------------------------------
LoadFullPath = sprintf('construction/RPTU_SCLDPC_Mb%d_Nb%d_Z%d.mat', mb, nb, Z);
profile = load(LoadFullPath);
BGCpd   = zeros(mb, nb, cpd_w); % construct basic matrices
BGCpd(:, :, 1) = profile.BGCpd(:, :, 1);
BGCpd(:, :, 2) = profile.BGCpd(:, :, 2);

% generate the coupled base graph for window
BGWin_len = winSize + cpd_w - 1;
BGWin     = zeros(mb*winSize, nb*BGWin_len) - 1;
for i_row = 1 : winSize
    for i_w = 1 : cpd_w
        BGWin((i_row-1)*mb+1 : i_row*mb, (i_w+i_row-1-1)*nb+1 : (i_w+i_row-1)*nb) = BGCpd(:, :, cpd_w-i_w+1);
    end
end

% generate an auxliary BG for termination 
BGTerm = BaseGraphTermAuxliaryMatrixGenerator(BGCpd, winSize, cpd_w, BGWin_len, mb, nb);

% generate an explicit H for each coupled component
HCpd = zeros(Mb, Nb, cpd_w);
for i_w = 1:cpd_w
    for ii = 1:mb
        for jj = 1:nb
            HCpd((ii-1)*Z+1:ii*Z, (jj-1)*Z+1:jj*Z, i_w) = identityMatrixExtend(Z, BGCpd(ii, jj, i_w));
        end
    end
end

% modify H0 as the standard format
[HCpd, IdentityTail, IdentityTailInv] = BinaryGaussianEliminationForSCLDPC(HCpd, Mb, Nb, cpd_w);

% extract H0'
H0Prime = HCpd(:, 1:Kb, 1);

% ---------------------------------------------------------------------------------
% encoding (parallel version)
% ---------------------------------------------------------------------------------
for i_run = 1 : frames : max_frames
    if nb_efer(end) >= max_errors && nb_ebler(end) >= 100
        break;
    end

    u = randi([0 1], Kb, frames, cpd_L);
    x = zeros(Nb, frames, cpd_L);
    partialSyndrome = zeros(Mb, frames, cpd_w-1); % partial syndrome registers: ps_1, ps_2, ps_3, ..., ps_w-1

    for i_L = 1 : cpd_L
        p = mod(H0Prime*u(:, :, i_L) + partialSyndrome(:, :, 1), 2) ;

        x(:, :, i_L) = [u(:, :, i_L); p];

        % implemenet cyclic shifting registers for partial syndromes
        partialSyndrome_last  = partialSyndrome;
        partailSyndrome_block = zeros(Mb, frames, cpd_w-1);
        for i_w = 1 : cpd_w - 1
            partailSyndrome_block(:, :, i_w) = mod(HCpd(:, :, i_w+1)*x(:, :, i_L), 2);
        end

        % cyclical shifting and updating
        for i_w = 1:cpd_w-2
            partialSyndrome(:, :, i_w) = mod(partialSyndrome_last(:, :, i_w+1) + partailSyndrome_block(:, :, i_w), 2);
        end
        partialSyndrome(:, :, end) = partailSyndrome_block(:, :, end);
    end

    % add a short tail composed of partial symdrome
    x = permute(x, [1,3,2]);

    x_sent_full = reshape(x, cpd_L*Nb, frames);
    for i_w = 1:cpd_w-1
        % post-processing for this tail (to avoid the influence of the above Gaussian elimination)
        partialSyndrome(:, :, i_w) = mod(IdentityTailInv*partialSyndrome(:, :, i_w), 2);

        x_sent_full = [x_sent_full; partialSyndrome(:, :, i_w)];
    end

    noise = randn(size(x_sent_full));
    for i_snr = 1:length(snr)
        if nb_efer(i_snr) >= max_errors && nb_ebler(i_snr) >= 100
            continue;
        end

        % ---------------------------------------------------------------------------------
        % channel (awgn, bpsk)
        % ---------------------------------------------------------------------------------
        sigma = 1/sqrt(2.0*rate) * 10^(-snr(i_snr)/20);

        symbol   = 1 - 2*x_sent_full;
        waveform = symbol + noise * sigma;
        chOut    = 2*waveform / sigma^2;
        
        % complement the input LLR vector to avoid the risk of array out of bounds
        chOut_full = [chOut; ones((winSize - 1)*Nb - (cpd_w-1)*Mb, frames)*1000];
        chOut_full_save = chOut_full;
        for i = 1:length(chOut_full_save(:))
            chOut_full_save(i) = (chOut_full_save(i) < 0);
        end 

        % reshape
        BGWin      = reshape(BGWin', length(BGWin(:)), 1);
        BGTerm     = reshape(BGTerm', length(BGTerm(:)), 1);
        if MultiCoreFlag == 0
            chOut_full = reshape(chOut_full, length(chOut_full(:)), 1);
            [rxcbs, rxitera] = sc_ldpc_layered_nms_float_decoding_mex(chOut_full, BGWin, BGTerm, winSize, cpd_L, cpd_w, Z, itera, nb, mb, Nb, Mb, alpha, frames);
            rxcbs = reshape(chOut_full, cpd_L*Nb, frames);
        else
            chOut_full = reshape(chOut_full, length(chOut_full(:))/groups, groups);
            rxcbs      = zeros(cpd_L*Nb*frames/groups, groups);
            rxitera    = zeros(1, groups);
            parfor i_groups = 1 : groups
                [rxcbs(:, i_groups), rxitera(i_groups)] = sc_ldpc_layered_nms_float_decoding_mex(chOut_full(:, i_groups), BGWin, BGTerm, winSize, cpd_L, cpd_w, Z, itera, nb, mb, Nb, Mb, alpha, frames/groups);
            end
            rxcbs   = reshape(rxcbs, cpd_L*Nb, frames);
            rxitera = sum(rxitera);
        end

        % compare the sent data with the decoded one
        dec_info_bits = zeros(Kb, cpd_L, frames);
        x_info_bits   = zeros(Kb, cpd_L, frames);
        x_sent_bits   = zeros(Kb, cpd_L, frames);
        for i_L = 1 : cpd_L-1
            dec_info_bits(:, i_L, :) = rxcbs((i_L-1)*Nb+1:(i_L-1)*Nb+Kb, :);
            x_info_bits(:, i_L, :)   = x_sent_full((i_L-1)*Nb+1:(i_L-1)*Nb+Kb, :);
            x_sent_bits(:, i_L, :)   = chOut_full_save((i_L-1)*Nb+1:(i_L-1)*Nb+Kb, :);
        end
        errNum  = sum(permute(abs(x_info_bits - dec_info_bits), [2 3 1]), 3);
        errSentNum  = sum(permute(abs(x_info_bits - x_sent_bits), [2 3 1]), 3);
        disp(['x_info_bits : ', num2str(sum(x_info_bits))]);
        disp(['x_sent_bits : ', num2str(sum(x_sent_bits))]);
        disp(['dec_info_bits : ', num2str(sum(dec_info_bits))]);
        disp(['errNum : ', num2str(sum(errNum))]);
        disp(['errSentNum : ', num2str(sum(errSentNum))]);

        % 
        nb_eber(i_snr)  = nb_eber(i_snr) + sum(errNum(:));
        nb_ebler(i_snr) = nb_ebler(i_snr) + sum(double(errNum(:) > 0));
        nb_efer(i_snr)  = nb_efer(i_snr) + sum(double(sum(errNum, 1) > 0));
        %nb_itera(i_snr) = nb_itera(i_snr) + rxitera;
        nb_frame(i_snr) = nb_frame(i_snr) + frames;
        
        nb_sent_ber(i_snr) = nb_sent_ber(i_snr) + sum(errSentNum(:));

        FER  = nb_efer./nb_frame;
        BLER = nb_ebler./nb_frame/cpd_L;
        BER  = nb_eber./nb_frame/cpd_L/Kb;
        ITER = nb_itera./nb_frame/cpd_L;

        %disp(' ');
        %disp(['nb_frame ' num2str(sum(nb_frame)) ' cpd_L ' num2str(cpd_L) ' Kb ' num2str(Kb)]);
        SENT = nb_sent_ber./nb_frame/cpd_L/Kb;
        %disp(['SENT ' num2str(SENT)]);

        framesNum = framesNum + frames;
        
        disp(' ');
        disp(['Running ' num2str(framesNum) ' frames']);
        %disp('Current block error performance');
        %disp(num2str([snr' FER' BLER' BER' ITER' nb_ebler' nb_efer']));
        toc
        break; % used to enhance the diversity of codewords
    end
end

%figure; semilogy(snr, FER, 'o-b', 'linewidth', 1); grid on;
%figure; semilogy(snr, BLER', 'o-b', 'linewidth', 1); grid on;
%figure; semilogy(snr, BER, 'o-r', 'linewidth', 1); grid on;
%xlabel('SNR (dB)');
%ylabel('Bit Error Rate (BER)');
%title('BER Performance');
%legend('Location', 'best');
%hold on


% Open a new figure window
figure;
semilogy(snr, SENT, 'o-r', 'LineWidth', 1);
hold on; % Hold on to plot additional data on the same axes
semilogy(snr, BER, 's-b', 'LineWidth', 1);
grid on;
xlabel('SNR (dB)');
ylabel('Bit Error Rate (BER)');
title('BER Performance Comparison');
set(gca, 'YScale', 'log');
ylim([1e-4 1]); 
xlim([2 6]);
legend('BER without decoding', 'BER with decoding', 'Location', 'best');
hold off;



%disp(['Dimensions of corrected_bits: ', num2str(size(decode_bits))]);
%disp(['corrected_bits: ', num2str(corrected_bits)]);
%Error using horzcat
%Dimensions of arrays being concatenated are not consistent.
%Error in spatially_coupled_ldpc_platform (line 237)
%disp(['corrected_bits: ', num2str(corrected_bits)]);

%disp(['Dimensions of errNum: ', num2str(size(errNum))]);
%disp(['Dimensions of correctedErrNum ', num2str(size(correctedErrNum))]);

%figure; semilogy(snr, CBER, 'o-r', 'linewidth', 1); grid on;
%xlabel('SNR (dB)');
%ylabel('Corrected Bit Error Rate (BER)');
%title('Corrected BER Performance');
%legend('Location', 'best');

%figure; semilogy(snr, correctedErrNum, 'o-r', 'linewidth', 1); grid on;
%xlabel('SNR (dB)');
%ylabel('Bit Error Rate (BER)');
%title('Corrected BER Performance');
%legend('Location', 'best');