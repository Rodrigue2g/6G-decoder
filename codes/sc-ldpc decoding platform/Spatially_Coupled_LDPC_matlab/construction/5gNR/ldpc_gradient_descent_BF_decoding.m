% ---------------------------------------------------------------------------------
% FUNCTION INFORMATION (c) 2023 Telecommunications Circuits Laboratory, EPFL
% ---------------------------------------------------------------------------------
% name  : ldpc_gradient_descent_BF_decoding
% descr : decoding of LDPC block codes (hard-decision) that support multiple modes

function [decode_bits, decode_itera] = ldpc_gradient_descent_BF_decoding(H, itera, LLRIn, frames)

% -------------------------------------------------------------------
% Parameter and memory allocations
% -------------------------------------------------------------------
[~, H_cols] = size(H);
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
    % Inversion/energy function
    % -------------------------------------------------------------------
    for i_itera = 1:itera
        for i_cols = 1:H_cols
            inverFunc(i_cols) = sum(full(H(:,i_cols)).*(1-2*syndrome)) + (1-2*xhat(i_cols))*LLRIn(i_cols, i_frames);
        end

        [~, flipIdx]  = min(inverFunc);
        xhat(flipIdx) = mod(xhat(flipIdx)+1, 2);
        syndrome      = mod(syndrome+H(:, flipIdx), 2);

        if sum(syndrome) == 0
            break;
        end
    end

    decode_bits(:, i_frames)  = xhat;
    decode_itera(:, i_frames) = i_itera;

end
end