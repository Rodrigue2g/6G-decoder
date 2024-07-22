% This function adds the 0 case
function [s] = signOP(x)
    s = 2*(x>=0) - 1;
end