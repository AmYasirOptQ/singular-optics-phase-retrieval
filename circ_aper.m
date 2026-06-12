function mask = circ_aper(N, R)
% CIRC_APER  Circular binary aperture mask.
%
%   mask = circ_aper(N, R)
%
%   N : grid size (pixels per side)
%   R : aperture radius (pixels)
%
%   Returns an N x N matrix: 1 inside the circle, 0 outside.
%   Circle is centred at (N/2+1, N/2+1) in 1-based indexing.

[cx, cy] = meshgrid(1:N, 1:N);
mask = double( (cx - N/2 - 1).^2 + (cy - N/2 - 1).^2 < R^2 );
end
