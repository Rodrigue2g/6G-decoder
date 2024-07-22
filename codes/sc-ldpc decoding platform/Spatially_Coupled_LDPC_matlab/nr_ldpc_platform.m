% -----------------------------------------------------------------------------------------------------------
% FUNCTION INFORMATION
% -----------------------------------------------------------------------------------------------------------
% name   : nr_ldpc_platform
% author : yuqing ren
% data   : 2023.09.18
% descr  : decoding 5g-nr ldpc with a series of algorithms (remove all special operations)

tic
rng(1);

% -----------------------------------------------------------------------------------------------------------
% 5G LDPC codes configuration
% -----------------------------------------------------------------------------------------------------------
TxRx.bgn = 1;   % base graph selection
TxRx.Z   = 192; 
TxRx.CR  = 1/2; % note that indiviual bit-flipping is incompatible with rate-matching (lack of necessary soft LLRs) 

if TxRx.bgn == 1
    TxRx.Nb = 68; % Fixed value for BG1
    TxRx.Mb = 46;
    TxRx.Kb = 22;
else
    TxRx.Nb = 52; % Fixed value FOR BG2
    TxRx.Mb = 42;
    TxRx.Kb = 10; % Directly set it as 10 to simplify the decoding
end

TxRx.usedMb  = ceil(TxRx.Kb/TxRx.CR) - TxRx.Kb;
TxRx.usedNb  = TxRx.usedMb + TxRx.Kb;

TxRx.N_graph = TxRx.Nb*TxRx.Z;
TxRx.M_graph = TxRx.Mb*TxRx.Z;
TxRx.K_graph = TxRx.Kb*TxRx.Z;

TxRx.K       = TxRx.K_graph;
TxRx.M       = TxRx.usedMb*TxRx.Z;

TxRx.padding = (TxRx.K + TxRx.M) - ceil(TxRx.K/TxRx.CR);
TxRx.CRTx    = TxRx.K/(TxRx.K+TxRx.M-TxRx.padding);

% -----------------------------------------------------------------------------------------------------------
% Modulation scheme settings
% -----------------------------------------------------------------------------------------------------------
TxRx.Qm = 2; % BPSK

% -----------------------------------------------------------------------------------------------------------
% Simulation settings
% -----------------------------------------------------------------------------------------------------------
TxRx.SNR_type  = 1; % 1 is Eb/N0
TxRx.SNRrange  = 1.0:0.1:2.0;

% Running setting
TxRx.group = 1;
TxRx.frame = 1000;
TxRx.disp  = 1000;
TxRx.max_iter = 10;                                 % Maximum number of iterations
TxRx.nb_frames_max = 1e9;                           % Maximum number of simulated frames
TxRx.nb_err_frames = 100*ones(size(TxRx.SNRrange)); % Maximum number of erroneous frames for each SNR value

% Extra optimizations
TxRx.beta = 1;
TxRx.norm = 3/4;

% -----------------------------------------------------------------------------------------------------------
% Launch simulation
% -----------------------------------------------------------------------------------------------------------
[rows, cols, ~, ~, ~] = getParameters(TxRx.bgn, TxRx.Z);
H = zeros(TxRx.M_graph, TxRx.N_graph); % Instantiate a H
for i_info = 1 : length(rows)
    H(rows(i_info), cols(i_info)) = 1;
end

% Note that here we only extract part of H based on the current M
% Since each entry of the prototype matrix is a ZxZ cyclic shift identity matrix or all zero matrix
% This pruning doesn't harm the sparsity of H, due to the independency of each block
H = H(1: TxRx.M, 1: TxRx.K_graph + TxRx.M); % For the current K and CR
ldpcencoder = comm.LDPCEncoder(sparse(H));  % Returns the position of all non-zero elements to decrease memory

% -----------------------------------------------------------------------------------------------------------
% Generation of the prototype matrix and parity-check matrix (base matrix of QC-LDPC)
% -----------------------------------------------------------------------------------------------------------
BG_2D     = get_BG(TxRx.bgn, TxRx.Z, H);
BG_r      = size(BG_2D,1);
BG_c      = size(BG_2D,2);
BG        = reshape(BG_2D', 1, BG_r*BG_c);
fileName  = sprintf('5gNR_BG%d_Mb%d_Nb%d_Z%d', TxRx.bgn, TxRx.usedMb, TxRx.usedNb, TxRx.Z);
full_path = fullfile('construction\ParityCheckMatrixProfiles', fileName);
save(full_path, 'BG_2D');

H_2D = zeros(BG_r*TxRx.Z, BG_c*TxRx.Z);
for ii = 1:BG_r
    for jj = 1:BG_c
        H_2D((ii-1)*TxRx.Z+1:ii*TxRx.Z, (jj-1)*TxRx.Z+1:jj*TxRx.Z) = identityMatrixExtend(TxRx.Z, BG_2D(ii, jj));
    end
end
H_2D = sparse(H_2D);

% -----------------------------------------------------------------------------------------------------------
% Prepare memory and variable
% -----------------------------------------------------------------------------------------------------------
SNR_list = TxRx.SNRrange;

% Instantiate the memory
nb_EF    = zeros(1, length(SNR_list));
nb_EB    = zeros(1, length(SNR_list));
nb_frame = zeros(1, length(SNR_list));
acc_iter = zeros(1, length(SNR_list));

FER = zeros(1, length(SNR_list));
BER = zeros(1, length(SNR_list));
average_iter = zeros(1, length(SNR_list));

x_full_seq = zeros(size(H,2), TxRx.frame);

% -----------------------------------------------------------------------------------------------------------
% Start the simulation
% -----------------------------------------------------------------------------------------------------------
for i_run = 1 : TxRx.frame : TxRx.nb_frames_max
    % Check the termination
    if nb_EF(end) >= TxRx.nb_err_frames
        disp(' ');
        disp(['Sim iteration running = ' num2str(nb_frame(end))]);
        disp(['BGN = ' num2str(TxRx.bgn) ' N = ' num2str(TxRx.K+TxRx.M-TxRx.padding) ' K = ' num2str(TxRx.K) ' CR = ' num2str(TxRx.CRTx) ' Z = ' num2str(TxRx.Z) ' Itera = ' num2str(TxRx.max_iter)]);
        disp(['EbN0 or SNR = ' num2str(TxRx.SNR_type)]);
        disp('Current block error performance');
        disp(num2str([SNR_list'  FER'  BER' average_iter' nb_EF']));
        disp(' ')
    end

    % Generate the source messages
    info_bits  = randi([0 1], TxRx.K, TxRx.frame);
    for i_frame = 1 : TxRx.frame
        x_full_seq(:, i_frame) = ldpcencoder(info_bits(:, i_frame)); % full codeword, length is (TxRx.K_graph + M)
    end

    x_sent_seq = x_full_seq;

    % Generate the noise
    noise = randn(size(x_sent_seq));

    % Check the termination
    if nb_EF(end) >= TxRx.nb_err_frames
        break;
    end

    for i_SNR = 1: length(SNR_list)
        if nb_EF(i_SNR) >= TxRx.nb_err_frames
            continue;
        end

        nb_frame(i_SNR) = nb_frame(i_SNR) + TxRx.frame;

        % Noise calculation comes from Prof. Wang's codes (in the folder called channel modulation)
        if TxRx.SNR_type == 1
            sigma = 1/sqrt(2.0*TxRx.CRTx*log2(TxRx.Qm)) * 10^(-SNR_list(i_SNR)/20);
        else
            sigma = 1/sqrt(2.0) * 10^(-SNR_list(i_SNR)/20);
        end

        symbols  = 1 - 2*x_sent_seq;
        waveform = symbols + noise * sigma;
        chOut    = 2*waveform / sigma^2;
        chOut(1:TxRx.padding, :) = 0; % rate-compatibility by puncturing

        % [rxcbs, actualniters] = ldpc_weighted_BF_decoding(H_2D, TxRx.max_iter, chOut*(sigma^2/2), TxRx.frame);
        % [rxcbs, actualniters] = ldpc_gradient_descent_BF_decoding(H_2D, TxRx.max_iter, chOut*(sigma^2/2), TxRx.frame);
        % [rxcbs, actualniters] = ldpc_multiple_gradient_descent_BF_decoding(H_2D, TxRx.max_iter, chOut*(sigma^2/2), TxRx.frame, -2);

        % Decode based on groups
        chOut = reshape(chOut, size(chOut, 1)*size(chOut, 2)/TxRx.group, TxRx.group);
        rxcbs = zeros(TxRx.group, size(chOut, 1));
        actualniters = zeros(TxRx.group, 1);
        parfor i_group = 1:TxRx.group
            [rxcbs(i_group, :), actualniters(i_group)] = ldpc_layered_nms_float_mex(BG, TxRx.norm, TxRx.Z, TxRx.max_iter, BG_r, BG_c, TxRx.frame/TxRx.group, chOut(:, i_group)', 1:BG_r);
        end
        rxcbs = reshape(rxcbs', 1, size(rxcbs,1)*size(rxcbs,2));
        actualniters = sum(actualniters);

        toc
        acc_iter(i_SNR) = acc_iter(i_SNR) + actualniters;
        rxcbs  = reshape(rxcbs, size(x_full_seq));
        rxcbsi = rxcbs(1: TxRx.K, :);

        % Compare the sent data with the decoded one
        for i_frame = 1: TxRx.frame
            if ~isequal(rxcbsi(:, i_frame), x_full_seq(1: TxRx.K, i_frame))
                nb_EF(i_SNR) = nb_EF(i_SNR) + 1;
                nb_EB(i_SNR) = nb_EB(i_SNR) + sum(rxcbsi(:, i_frame) ~= x_full_seq(1: TxRx.K, i_frame));
            end
        end

        % Calculate the FER/BER/Itera
        FER = nb_EF./nb_frame;
        BER = nb_EB./(TxRx.K*nb_frame);
        average_iter = acc_iter./nb_frame;

        % skip
        break;
    end

    % disp
    if  mod(i_run + TxRx.frame, TxRx.disp) == 1
        disp(' ');
        disp(['Sim iteration running = ' num2str(sum(nb_frame))]);
        disp(['BGN = ' num2str(TxRx.bgn) ' N = ' num2str(TxRx.K+TxRx.M-TxRx.padding) ' K = ' num2str(TxRx.K) ' CR = ' num2str(TxRx.CRTx) ' Z = ' num2str(TxRx.Z) ' Itera = ' num2str(TxRx.max_iter)]);
        disp(['EbN0 or SNR = ' num2str(TxRx.SNR_type)]);
        disp('Current block error performance');
        disp(num2str([SNR_list'  FER'  BER' average_iter' nb_EF']));
        disp(' ')
    end

    if nb_EF(end) >= TxRx.nb_err_frames
        disp(' ');
        disp(['Sim iteration running = ' num2str(sum(nb_frame))]);
        disp(['BGN = ' num2str(TxRx.bgn) ' N = ' num2str(TxRx.K+TxRx.M-TxRx.padding) ' K = ' num2str(TxRx.K) ' CR = ' num2str(TxRx.CRTx) ' Z = ' num2str(TxRx.Z) ' Itera = ' num2str(TxRx.max_iter)]);
        disp(['EbN0 or SNR = ' num2str(TxRx.SNR_type)]);
        disp('Current block error performance');
        disp(num2str([SNR_list'  FER'  BER' average_iter' nb_EF']));
        disp(' ')
        break;
    end
end

figure; semilogy(SNR_list, FER, 'o-b', 'linewidth', 1); grid on;
figure; semilogy(SNR_list, BER, 'o-r', 'linewidth', 1); grid on;
