% -----------------------------------------------------------------------------------------------------------
% FUNCTION INFORMATION
% -----------------------------------------------------------------------------------------------------------
% name   : ldpc_layered_ms_float
% author : yuqing ren
% data   : 2022.05.28
% description : layered MS decoding algorithm for the floating performance (block parallel, log domain)

function [decode_bits, decode_itera] = ldpc_layered_ms_float(BG, Z, itera, LLRIn, layer_order)

% -------------------------------------------------------------------
% Parameter and memory allocations
% -------------------------------------------------------------------
[H_rows, H_cols] = size(BG);

QMem = reshape(LLRIn, Z, H_cols); % Q memory: store the soft value of each VN
RMem = zeros(Z, H_rows*H_cols);   % R memory: store the internal message from CN to VN

% Temporary memory
TMem    = zeros(Z, H_cols);       % T memory: store the internal temporary message during the layered calculation
QSign   = zeros(Z, H_cols);       % Record the sign bit of Q memory

% -------------------------------------------------------------------
% Decoding
% -------------------------------------------------------------------
for i_itera = 1 : itera
    TDFlag   = 1; % syndrome check
    
    for i_rows = 1 : H_rows
        SynCheck = zeros(Z, 1);
        
% -------------------------------------------------------------------
% Phase 1: MIN
% -------------------------------------------------------------------
        min1   = Inf * ones(Z, 1); % store the first minima
        min2   = Inf * ones(Z, 1); % store the second minima
        sign   = ones(Z, 1);       % store the sign bit of TMem (Z-based)

        minIdx = zeros(Z, 1); % store the column index corresponding to min1
        minVec = zeros(Z, 1); % store the minima except the current column

        for i_cols = 1 : H_cols
            % checkout whether there is -1 in base matrix, if yes, directly skip
            if BG(layer_order(i_rows), i_cols) == -1
                continue;
            end
            shiftNum = BG(layer_order(i_rows), i_cols); % cyclic shifting value (c of P^c)
            
            % Calculate the Qc and TMem
            Qc = Rotation(QMem(:, i_cols), shiftNum, Z); % rotate left by shiftNum bits
            TMem(:, i_cols) = Qc - RMem(:, (layer_order(i_rows)-1)*H_cols + i_cols); % TMem is actually the previous Qvc (from VN --> CN)
            
            % Get the sign and magnitude
            sign             = sign.*signOP(TMem(:, i_cols));
            TMemAbs          = abs(TMem(:, i_cols));
            QSign(:, i_cols) = signOP(Qc);
            
            % Search the min1, min2, and the column index corresponding to min1
            for i_Z = 1 : Z
                if TMemAbs(i_Z) < min1(i_Z)
                    min2(i_Z)   = min1(i_Z);
                    min1(i_Z)   = TMemAbs(i_Z);
                    minIdx(i_Z) =  i_cols;
                elseif TMemAbs(i_Z) < min2(i_Z)
                    min2(i_Z)   = TMemAbs(i_Z);
                end
            end
        end
        
% -------------------------------------------------------------------
% Phase 2: Q- and R- messages update
% -------------------------------------------------------------------
        for i_cols = 1 : H_cols
            % checkout whether there is -1 in base matrix, if yes, directly skip
            if BG(layer_order(i_rows), i_cols) == -1
                continue;
            end
            shiftNum = BG(layer_order(i_rows), i_cols); % cyclic shifting value (c of P^c)
            
            % Deal with the special cases that the minina belongs to the current column
            for i_Z = 1 : Z
                if minIdx(i_Z) == i_cols
                    minVec(i_Z) = min2(i_Z);
                else
                    minVec(i_Z) = min1(i_Z);
                end
            end
            
            RMem(:, (layer_order(i_rows)-1)*H_cols + i_cols) = sign.*signOP(TMem(:, i_cols)).*minVec;
            Qtmp = TMem(:, i_cols) + RMem(:, (layer_order(i_rows)-1)*H_cols + i_cols);
            
% -------------------------------------------------------------------
% Phase 3: Syndrome check
% -------------------------------------------------------------------
            % Syndrome check (execute the verification by each entry in BG)
            if sum(abs((Qtmp < 0) - (QSign(:, i_cols) < 0))) ~= 0
                TDFlag = 0;
            end
            SynCheck = xor(SynCheck, (Qtmp < 0));
            
            % Inverse rotation, i.e., P^c*P^c' = I
            Qtmp = Rotation(Qtmp, Z - shiftNum, Z);
            QMem(:, i_cols) = Qtmp;
        end
        
        if sum(SynCheck) > 0
            TDFlag = 0;
        end
    end
    
    % check whether TDFlag == 1 for each iteration
    if TDFlag == 1
        break;
    end
end

% generate the output
decode_bits  = double(reshape(QMem, H_cols*Z, 1) < 0);
decode_itera = i_itera; 
end

% This function could perform the cyclic shifting
function yOut = Rotation(yIn, shiftNum, Z)
    yOut = [yIn(shiftNum + 1 : Z); yIn(1 : shiftNum)];
end
