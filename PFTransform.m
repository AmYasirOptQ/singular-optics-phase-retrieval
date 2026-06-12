function u_out = PFTransform(u_in, N)
% PFTRANSFORM  Forward partial Fourier transform.
%
%   u_out = PFTransform(u_in, N)
%
%   Maps a Laguerre-Gaussian field to a rotated Hermite-Gaussian field.
%   Does not conserve longitudinal OAM — can create or annihilate vortices.
%
%   Implementation: u_out = B2 * u_in
%   where B2 = fftshift(DFT_N) / sqrt(N)

[c, r] = meshgrid(0:N-1, 0:N-1);
B2     = fftshift(exp(-2i*pi/N) .^ (r.*c)) / sqrt(N);
u_out  = B2 * u_in;
end
