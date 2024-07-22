% ---------------------------------------------------------------------------------
% FUNCTION INFORMATION (c) 2023 Telecommunications Circuits Laboratory, EPFL
% ---------------------------------------------------------------------------------
% name  : ldpc_weighted_BF_decoding
% descr : decoding of LDPC block codes (hard-decision)

function [decode_bits, decode_itera] = ldpc_weighted_BF_decoding(H, itera, LLRIn, frames)

% -------------------------------------------------------------------
% Parameter and memory allocations
% -------------------------------------------------------------------
[H_rows, H_cols] = size(H);
weight    = zeros(1, H_rows);
inverFunc = zeros(1, H_cols);

decode_bits  = zeros(size(LLRIn, 1), frames);
decode_itera = zeros(1, frames);

for i_frames = 1:frames
    % -------------------------------------------------------------------
    % Synchrome Calculation
    % -------------------------------------------------------------------
    xhat     = double(LLRIn(:, i_frames) < 0); % hard decision
    syndrome = mod(H*xhat, 2); % synchrome calculation

    % -------------------------------------------------------------------
    % Compute the weights
    % -------------------------------------------------------------------
    for i_rows = 1:H_rows
        vectorTmp      = full(H(i_rows, :)'.*abs(LLRIn(:, i_frames)));
        nonZeroTmp     = vectorTmp(vectorTmp ~= 0);
        weight(i_rows) = min(nonZeroTmp);
    end

    % -------------------------------------------------------------------
    % Inversion function
    % -------------------------------------------------------------------
    for i_itera = 1:itera
        for i_cols = 1:H_cols
            inverFunc(i_cols) = sum(full(H(:,i_cols)).*weight'.*(2*syndrome-1));
        end

        [~, flipIdx]  = max(inverFunc);
        xhat(flipIdx) = mod(xhat(flipIdx)+1, 2);
        syndrome      = mod(syndrome+H(:, flipIdx),2);

        if sum(syndrome) == 0
            break;
        end
    end

    decode_bits(:, i_frames)  = xhat;
    decode_itera(:, i_frames) = i_itera;

end
end