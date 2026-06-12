function L = laguerre(p, alpha, x)
% LAGUERRE  Generalised Laguerre polynomial L_p^alpha(x).
% Three-term recurrence relation; x may be a matrix.
if p == 0
    L = ones(size(x));
    return
end
L_prev = ones(size(x));
L_curr = 1 + alpha - x;
for n = 2 : p
    L_next = ((2*n - 1 + alpha - x) .* L_curr ...
             - (n - 1 + alpha)      .* L_prev) / n;
    L_prev = L_curr;
    L_curr = L_next;
end
L = L_curr;
end
