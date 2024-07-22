% ---------------------------------------------------------------------------------
% FUNCTION INFORMATION (c) 2023 Telecommunications Circuits Laboratory, EPFL
% ---------------------------------------------------------------------------------
% name  : identityMatrixExtend
% descr : transform BG to H
 
function qout = identityMatrixExtend(Z, cyclicValue)
    if cyclicValue == -1
        qout = zeros(Z);
    else
        qout = circshift(eye(Z), [0 cyclicValue]);
    end
end