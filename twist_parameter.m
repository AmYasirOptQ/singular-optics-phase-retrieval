function tau = twist_parameter(u, X, Y, dx, lambda)
% TWIST_PARAMETER  Orbital angular momentum twist parameter.
%
%   tau = twist_parameter(u, X, Y, dx, lambda)
%
%   Computes  tau = (<x py> - <y px>) / hbar
%                 - (<x><py> - <y><px>) / hbar
%
%   where hbar = lambda/(2*pi) and momenta use sixth-order finite differences.
%   For a pure LG_p^l mode, tau = l (the topological charge).

lambdabar = lambda / (2*pi);

[dudx, dudy] = third_order_gradient(u, dx);

uc = conj(u);
D  = sum(sum(uc .* u));

x_m  = sum(sum(uc .* X    .* u))   / D;
y_m  = sum(sum(uc .* Y    .* u))   / D;
px_m = sum(sum(uc .* dudx)) * (-1i*lambdabar) / D;
py_m = sum(sum(uc .* dudy)) * (-1i*lambdabar) / D;

xpy  = sum(sum(uc .* X .* dudy)) * (-1i*lambdabar) / D;
ypx  = sum(sum(uc .* Y .* dudx)) * (-1i*lambdabar) / D;

tau  = (xpy - ypx) / lambdabar ...
     - (x_m*py_m - y_m*px_m) / lambdabar;
end
