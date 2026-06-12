function u_out = fresnel_propagation(u_in, z, lambda, L)
% FRESNEL_PROPAGATION  Free-space propagation by the angular-spectrum method.
%
%   u_out = fresnel_propagation(u_in, z, lambda, L)
%
%   u_in   : N x N complex input field
%   z      : propagation distance (m); negative value = backward propagation
%   lambda : wavelength (m)
%   L      : screen size (m)
%
%   Transfer function:  H = exp(-i*pi*lambda*z*(fx^2+fy^2))

N  = size(u_in, 1);
dx = L / N;
df = 1 / (N*dx);

f  = (-N/2 : N/2-1) * df;
[fx, fy] = meshgrid(f, f);

H     = exp(-1i*pi*lambda*z*(fx.^2 + fy.^2));
U     = fftshift(fft2(fftshift(u_in))) * dx^2;
u_out = ifftshift(ifft2(ifftshift(U.*H))) * (N*df)^2;
end
