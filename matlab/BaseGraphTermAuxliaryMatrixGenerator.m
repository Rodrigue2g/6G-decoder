% ---------------------------------------------------------------------------------
% FUNCTION INFORMATION (c) 2023 Telecommunications Circuits Laboratory, EPFL
% ---------------------------------------------------------------------------------
% name  : BGAuxilaryGeneration
% descr : generate the auxiliary basegraph for the last terminated bits

function [BaseGraphTermMatrix] = BaseGraphTermAuxliaryMatrixGenerator(BaseGraph, winSize, cpd_w, BaseGraphWin_len, mb, nb)
% parameter definition
winDelay = winSize - 1; % delay of sliding in and sliding out

BaseGraphTerm_cols  = BaseGraphWin_len + winDelay - 1;
BaseGraphTerm_rows  = winSize + winDelay - 1;
BaseGraphTermMatrix = zeros(BaseGraphTerm_rows*mb, BaseGraphTerm_cols*nb) - 1;

% fill the component matrices
for i_rows = 1 : BaseGraphTerm_rows
    for i_w = 1 : cpd_w
        BaseGraphTermMatrix((i_rows-1)*mb+1 : i_rows*mb, (i_rows+i_w-2)*nb+1 : (i_rows+i_w-1)*nb) = BaseGraph(:, :, cpd_w-i_w+1);
    end
end
BaseGraphTermMatrix(:, (BaseGraphWin_len-1)*nb+1 : end) = -1;

% add an identity matrix
IdentityTail = eye(mb*(cpd_w-1)) - 1;
BaseGraphTermMatrix((winSize-1)*mb+1 : (winSize-1)*mb+mb*(cpd_w-1), (BaseGraphWin_len-1)*nb+1 : (BaseGraphWin_len-1)*nb+(cpd_w-1)*mb) = IdentityTail;

end

