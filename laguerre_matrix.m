function L = laguerre_matrix(p, alpha, x)
% LAGUERRE_MATRIX  Element-wise generalised Laguerre polynomial on a 2-D array.
% Wrapper around laguerre() for naming consistency with the paper.
L = laguerre(p, alpha, x);
end
