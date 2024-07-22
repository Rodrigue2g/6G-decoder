% ---------------------------------------------------------------------------------
% FUNCTION INFORMATION (c) 2023 Telecommunications Circuits Laboratory, EPFL
% ---------------------------------------------------------------------------------
% name  : ldpc_multiple_gradient_descent_BF_decoding
% descr : decoding of LDPC block codes (hard-decision) that support multiple modes

function [decode_bits, decode_itera] = ldpc_multiple_gradient_descent_BF_decoding(H, itera, LLRIn, frames, delta)

% -------------------------------------------------------------------
% Parameter and memory allocations
% -------------------------------------------------------------------
[~, H_cols] = size(H);
inverFunc = zeros(1, H_cols);

decode_bits  = zeros(size(LLRIn, 1), frames);
decode_itera = zeros(1, frames);

for i_frames = 1:frames
    % -------------------------------------------------------------------
    % Synchrome/Objective function Calculation
    % -------------------------------------------------------------------
    xhat     = double(LLRIn(:, i_frames) < 0); % hard decision
    syndrome = mod(H*xhat, 2); % synchrome calculation

    objectFunc = (1-2*xhat)'*LLRIn(:, i_frames) + sum((1-2*syndrome));
    modeFlag   = 1;

    % -------------------------------------------------------------------
    % Inversion/energy function
    % -------------------------------------------------------------------
    for i_itera = 1:itera
        for i_cols = 1:H_cols
            inverFunc(i_cols) = sum(full(H(:,i_cols)).*(1-2*syndrome)) + (1-2*xhat(i_cols))*LLRIn(i_cols, i_frames);
        end

        [values, Idx]  = sort(inverFunc, 'ascend');
        
        if modeFlag == 1
            flipIdx = Idx(values < delta);
            if isempty(flipIdx)
                flipIdx   = Idx(1);
                modeFlag  = 0;
            end
            xhat(flipIdx) = mod(xhat(flipIdx)+1, 2);
            syndrome      = mod(syndrome+sum(H(:, flipIdx),2), 2);

        else
            flipIdx       = Idx(1);
            xhat(flipIdx) = mod(xhat(flipIdx)+1, 2);
            syndrome      = mod(syndrome+H(:, flipIdx), 2);
        end

        % Update objective function
        objectLast = objectFunc;
        objectFunc = (1-2*xhat)'*LLRIn(:, i_frames) + sum((1-2*syndrome));

        if objectFunc < objectLast
            modeFlag = 0;
        end

        if sum(syndrome) == 0
            break;
        end
    end

    decode_bits(:, i_frames)  = xhat;
    decode_itera(:, i_frames) = i_itera;

end
end