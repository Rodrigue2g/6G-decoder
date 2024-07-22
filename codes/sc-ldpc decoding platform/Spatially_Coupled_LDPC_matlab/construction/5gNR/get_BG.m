% -----------------------------------------------------------------------------------------------------------
% FUNCTION INFORMATION
% -----------------------------------------------------------------------------------------------------------
% name   : gen_Hp
% author : yuqing ren
% data   : 2022.05.28
% descr  : generate the correct BG (basic graph)

function BG = get_BG(bgn, Z, H)

% Pre-generated 5G lifting sizes
lifting_sizes = [2,  4,  8,  16,  32,  64,  128, 256, ...
                 3,  6,  12, 24,  48,  96,  192, 384, ...
                 5,  10, 20, 40,  80,  160, 320, ...
                 7,  14, 28, 56,  112, 224, ...
                 9,  18, 36, 72,  144, 288, ...
                 11, 22, 44, 88,  176, 352, ...
                 13, 26, 52, 104, 208, ...
                 15, 30, 60, 120, 240];
set_index = [1, 1, 1, 1, 1, 1, 1, 1, ...
             2, 2, 2, 2, 2, 2, 2, 2, ...
             3, 3, 3, 3, 3, 3, 3, ...
             4, 4, 4, 4, 4, 4, ...
             5, 5, 5, 5, 5, 5, ...
             6, 6, 6, 6, 6, 6, ...
             7, 7, 7, 7, 7, ...
             8, 8, 8, 8, 8];
index = set_index(find(lifting_sizes == Z));

% load the BG matrix
if bgn == 1
    load parity_check_matrices_protocol_1
    BG = parity_check_matrices_protocol_1(:, :, index);
elseif bgn == 2
    load parity_check_matrices_protocol_2
    BG = parity_check_matrices_protocol_2(:, :, index);
else
    error('wrong base graph index in ldpc encoding.');
end

% Modulo with the current Z
for i_row = 1: length(BG(:, 1))
    for i_column = 1: length(BG(1, :))
        if (BG(i_row, i_column) > -1)
            BG(i_row, i_column) = mod(BG(i_row, i_column), Z); % mod operation for different Z
        end
    end
end

% Perform the pruning as same as H
BG = BG(1:length(H(:,1))/Z, 1:length(H(1,:))/Z);
end