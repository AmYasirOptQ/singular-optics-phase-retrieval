function u_out = IPFTransform(u_in, N)
% IPFTRANSFORM  Inverse partial Fourier transform.
%
%   u_out = IPFTransform(u_in, N)
%
%   Inverts PFTransform via conjugate transpose:  u_out = B2' * u_in
%   NaN values from numerical noise are zeroed.

[c, r] = meshgrid(0:N-1, 0:N-1);
B2     = fftshift(exp(-2i*pi/N) .^ (r.*c)) / sqrt(N);
u_out  = B2' * u_in;
u_out(isnan(u_out)) = 0;
end
