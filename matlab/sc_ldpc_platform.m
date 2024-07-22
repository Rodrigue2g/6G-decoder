% ---------------------------------------------------------------------------------
% FUNCTION INFORMATION (c) 2024 Telecommunications Circuits Laboratory, EPFL
% ---------------------------------------------------------------------------------
% name  : MultiFrame_TB
% descr : Testbench used to characterize the performance of bit flipping algorithms
%         by plotting BER curves.  

% ---------------------------------------------------------------------------------
% parameter definition
% ---------------------------------------------------------------------------------

tic
rng(1);
cpd_L   = 50;
cpd_w   = 2; % also defined as ms or mcc = w-1

nb  = 40;
mb  = 8;
Z   = 16;
kb  = nb-mb;

Nb = nb*Z;
Mb = mb*Z;
Kb = kb*Z;


snr   = 4:0.25:11; %2:0.25:8.5;    5.3 3.1 4 4.5 5.7   % change this to iterate over desired SNR values
frames = 100;             % change this to iterate over multiple frames (keep at 1 for now else it too damn slow)
nb_algorithms = 10; %7     % 2 (1 for the baseline and 1 for the algorithm to test) change this if we want to test more algorithms (see bellow)

% (cpd_w-1)*Mb is the length of tail we attach an identity matrix 
% at the end to do the termination
rate  = cpd_L*Kb / (cpd_L*Nb + (cpd_w-1)*Mb);


% list that will store all nb of errors
errors=zeros(nb_algorithms,length(snr), cpd_L, frames);
rxcbs =zeros(nb_algorithms, cpd_L*Nb);


% ---------------------------------------------------------------------------------
% Hw generation
% ---------------------------------------------------------------------------------

LoadFullPath = sprintf('construction/RPTU_SCLDPC_Mb8_Nb40_Z16', mb, nb, Z);
profile = load(LoadFullPath);
BGCpd   = zeros(mb, nb, cpd_w); % construct basic matrices
BGCpd(:, :, 1) = profile.BGCpd(:, :, 1);
BGCpd(:, :, 2) = profile.BGCpd(:, :, 2);


%%--------------- UNCOMMENT IF WE WANT TO TEST THE WINDOWED DECODER -------------
% CreateBG_WD_decoder;
   
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
[HCpd, ~, IdentityTailInv] = BinaryGaussianEliminationForSCLDPC(HCpd, Mb, Nb, cpd_w);

% extract H0'
H0Prime = HCpd(:, 1:Kb, 1);

% Generate the full H matrix
BaseGraphFull = zeros(mb*(cpd_L+cpd_w-1), nb*cpd_L + (cpd_w-1)*mb) - 1;
for i_cols = 1 : cpd_L
    for i_w = 1 : cpd_w
        BaseGraphFull((i_cols+i_w-2)*mb+1 : (i_cols+i_w-1)*mb, (i_cols-1)*nb+1 : i_cols*nb) = BGCpd(:, :, i_w);
    end
end
IdentityTail = eye(mb*(cpd_w-1)) - 1;
BaseGraphFull(cpd_L*mb+1 : end, cpd_L*nb+1 : end) = IdentityTail;

HFull = zeros(size(BaseGraphFull,1)*Z, size(BaseGraphFull,2)*Z);
for ii = 1 : size(BaseGraphFull,1)
    for jj = 1 : size(BaseGraphFull,2)
        HFull((ii-1)*Z+1:ii*Z, (jj-1)*Z+1:jj*Z) = identityMatrixExtend(Z, BaseGraphFull(ii, jj));
    end
end

for i_frame = 1:frames

    % random messages
    u = randi([0 1], Kb, cpd_L);

    % ---------------------------------------------------------------------------------
    % encoding (parallel version)
    % ---------------------------------------------------------------------------------
    x = zeros(Nb, cpd_L);
    partialSyndrome = zeros(Mb, cpd_w-1); % partial syndrome registers: ps_1, ps_2, ps_3, ..., ps_w-1

    for i_L = 1 : cpd_L
        p = mod(H0Prime*u(:, i_L) + partialSyndrome(:, 1), 2) ;

        x(:, i_L) = [u(:, i_L); p];

        % implement cyclic shifting registers for partial syndromes
        partialSyndrome_last  = partialSyndrome;
        partailSyndrome_block = zeros(Mb, cpd_w-1);
        for i_w = 1 : cpd_w - 1
            partailSyndrome_block(:, i_w) = mod(HCpd(:, :, i_w+1)*x(:, i_L), 2);
        end 

        % cyclical shifting and updating
        for i_w = 1:cpd_w-2
            partialSyndrome(:, i_w) = mod(partialSyndrome_last(:, i_w+1) + partailSyndrome_block(:, i_w), 2);
        end
        partialSyndrome(:, end) = partailSyndrome_block(:, end);
    end

    % add a short tail composed of partial symdrome
    x_sent_full = reshape(x, [numel(x),1]);
    for i_w = 1:cpd_w-1
        % post-processing for this tail (to avoid the influence of the above Gaussian elimination)
        partialSyndrome(:, i_w) = mod(IdentityTailInv*partialSyndrome(:, i_w), 2);
        x_sent_full = [x_sent_full; partialSyndrome(:, i_w)];
    end

    % the vector to compare to the decoded one
    x_info_bits   = zeros(Kb, cpd_L);
    for i_L = 1 : cpd_L-1
        x_info_bits(:, i_L)   = x_sent_full((i_L-1)*Nb+1:(i_L-1)*Nb+Kb, :);
    end    
    
    % keep the same random noise through all the SNR's
    noise = randn(size(x_sent_full));

    % iterate over all SNRs
    num_flips_all = zeros(nb_algorithms, length(snr), 100);  % Pre-allocate array for number of flips
    for i_snr = 1:length(snr)

        % ---------------------------------------------------------------------------------
        % channel (awgn, bpsk)
        % ---------------------------------------------------------------------------------

        sigma = 1/sqrt(2.0*rate) * 10^(-snr(i_snr)/20);
        chOut = (1 - 2*x_sent_full)+noise*sigma;
        chOut_neg = (2*x_sent_full - 1) + noise*sigma;
    
        % ---------------------------------------------------------------------------------
        % decoding 
        % ---------------------------------------------------------------------------------
        rxcbs(1,:)  = baseline_mex(HFull', chOut, cpd_L, cpd_w, Z, nb, mb);
        %rxcbs(2,:)  = SGDBF_mex(HFull', chOut, cpd_L, cpd_w, Z, nb, mb);   
        rxcbs(2,:)  = SGDBF_V2_0_mex(HFull', chOut, cpd_L, cpd_w, Z, nb, mb);
        rxcbs(3,:)  = SGDBF_V2_1_mex(HFull', chOut, cpd_L, cpd_w, Z, nb, mb);
        rxcbs(4,:)  = SGDBF_V2_2_mex(HFull', chOut, cpd_L, cpd_w, Z, nb, mb);
        rxcbs(5,:)  = SGDBF_V2_3_mex(HFull', chOut, cpd_L, cpd_w, Z, nb, mb);
        rxcbs(6,:)  = SGDBF_V2_4_mex(HFull', chOut, cpd_L, cpd_w, Z, nb, mb);
        rxcbs(7,:)  = SGDBF_V2_5_mex(HFull', chOut, cpd_L, cpd_w, Z, nb, mb);
        rxcbs(8,:)  = SGDBF_V2_6_mex(HFull', chOut, cpd_L, cpd_w, Z, nb, mb);
        rxcbs(9,:)  = SGDBF_V2_7_mex(HFull', chOut, cpd_L, cpd_w, Z, nb, mb);
        rxcbs(10,:)  = SGDBF_V2_8_mex(HFull', chOut, cpd_L, cpd_w, Z, nb, mb);

        %rxcbs(4,:)  = LC_HWBF_mex(HFull', chOut_neg, cpd_L, cpd_w, Z, nb, mb);

        %[rxcbs(1,:), num_flips] = baseline_mex(HFull', chOut, cpd_L, cpd_w, Z, nb, mb);
        %num_flips_all(1, i_snr, :) = num_flips;

        %[rxcbs(1,:), num_flips] = LC_HWBF_mex(HFull', chOut_neg, cpd_L, cpd_w, Z, nb, mb);
        %num_flips_all(1, i_snr, :) = num_flips;

        %[rxcbs(1,:), num_flips] = SGDBF_mex(HFull', chOut, cpd_L, cpd_w, Z, nb, mb);
        %num_flips_all(1, i_snr, :) = num_flips;

        %[rxcbs(2,:), num_flips] = SGDBF_V2_mex(HFull', chOut, cpd_L, cpd_w, Z, nb, mb);
        %num_flips_all(2, i_snr, :) = num_flips;
    
        %[rxcbs(4,:), num_flips] = MGDBF_mex(HFull', chOut, cpd_L, cpd_w, Z, nb, mb);
        %num_flips_all(4, i_snr, :) = num_flips;

        %rxcbs(2,:), num_flips] = MGDBF_V2_mex(HFull', chOut, cpd_L, cpd_w, Z, nb, mb);
        %num_flips_all(2, i_snr, :) = num_flips;



        %rxcbs(1,:) = baseline_mex(HFull', chOut, cpd_L, cpd_w, Z, nb, mb);
%
        %[rxcbs(2,:), ~] = LC_HWBF_mex(HFull', chOut_neg, cpd_L, cpd_w, Z, nb, mb);
%
        %[rxcbs(3,:), ~] = SGDBF_mex(HFull', chOut, cpd_L, cpd_w, Z, nb, mb);
        %[rxcbs(4,:), ~] = SGDBF_V2_mex(HFull', chOut, cpd_L, cpd_w, Z, nb, mb);
%
        %[rxcbs(5,:), ~] = MGDBF_mex(HFull', chOut, cpd_L, cpd_w, Z, nb, mb);
%
        %rxcbs(6,:)  = MGDBF_V2_1_mex(HFull', chOut, cpd_L, cpd_w, Z, nb, mb);
        %rxcbs(7,:)  = MGDBF_V2_2_mex(HFull', chOut, cpd_L, cpd_w, Z, nb, mb);
        %rxcbs(8,:)  = MGDBF_V2_3_mex(HFull', chOut, cpd_L, cpd_w, Z, nb, mb);

        %rxcbs(1,:)  = MGDBF_V2_1_mex(HFull', chOut, cpd_L, cpd_w, Z, nb, mb);
        %rxcbs(2,:)  = MGDBF_V2_2_mex(HFull', chOut, cpd_L, cpd_w, Z, nb, mb);
        %rxcbs(3,:)  = MGDBF_V2_3_mex(HFull', chOut, cpd_L, cpd_w, Z, nb, mb);
        %rxcbs(4,:)  = MGDBF_V2_4_mex(HFull', chOut, cpd_L, cpd_w, Z, nb, mb);
        %rxcbs(5,:)  = MGDBF_V2_5_mex(HFull', chOut, cpd_L, cpd_w, Z, nb, mb);
        %rxcbs(6,:)  = MGDBF_V2_6_mex(HFull', chOut, cpd_L, cpd_w, Z, nb, mb);
        %rxcbs(7,:)  = MGDBF_V2_7_mex(HFull', chOut, cpd_L, cpd_w, Z, nb, mb);
        %rxcbs(8,:)  = MGDBF_V2_8_mex(HFull', chOut, cpd_L, cpd_w, Z, nb, mb);
        %rxcbs(9,:)  = MGDBF_V2_9_mex(HFull', chOut, cpd_L, cpd_w, Z, nb, mb);
       % rxcbs(6,:)  = MGDBF_V3_mex(HFull', chOut, cpd_L, cpd_w, Z, nb, mb);
       

        % ---------------------------------------------------------------------------------
        % Find Nb of errors
        % ---------------------------------------------------------------------------------
        for code = 1:nb_algorithms
            % extract info bits
            dec_info_bits = zeros(Kb, cpd_L);
            for i_L = 1 : cpd_L-1
                dec_info_bits(:, i_L) = rxcbs(code, (i_L-1)*Nb+1:(i_L-1)*Nb+Kb);
            end

            % find the number of errors
            for i_L = 1:cpd_L
                errors(code, i_snr, i_L, i_frame)=sum(abs(x_info_bits(:,i_L)-dec_info_bits(:,i_L)),"all");
            end
        end

        disp(['Frame ' num2str(i_frame) ' at SNR ' num2str(snr(i_snr)) ' dB is done'])

    end
end

toc

% ---------------------------------------------------------------------------------
% Results
% ---------------------------------------------------------------------------------

view_data