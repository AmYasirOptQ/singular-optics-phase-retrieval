function [X, Y, r, phi, dx] = make_grid(N, L)
% MAKE_GRID  Spatial grid for paraxial beam simulation.
%   [X, Y, r, phi, dx] = make_grid(N, L)
%   N  : number of pixels per side
%   L  : total screen size (m)
dx  = L / N;
x   = (-N/2 : N/2-1) * dx;
[X, Y] = meshgrid(x, x);
r   = sqrt(X.^2 + Y.^2);
phi = atan2(Y, X);
end
