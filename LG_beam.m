function u = LG_beam(r, phi, l, p, w1, z1, zr, R1, lambda)
% LG_BEAM  Normalised Laguerre-Gaussian mode LG_p^l including propagation phase.
%
%   u = LG_beam(r, phi, l, p, w1, z1, zr, R1, lambda)
%
%   r, phi  : polar coordinate matrices (from make_grid or offset grid)
%   l, p    : azimuthal and radial indices
%   w1      : beam radius at evaluation plane  w1 = w0*sqrt(1+(z1/zR)^2)
%   z1      : propagation distance from waist (m)
%   zr      : Rayleigh range  zR = pi*w0^2/lambda
%   R1      : wavefront radius  R1 = z1*(1+(zR/z1)^2)
%   lambda  : wavelength (m)

k = 2*pi / lambda;

% amplitude envelope
A = (1/w1) * sqrt( (4*factorial(p)) / (pi*2*factorial(abs(l)+p)) );
u = A .* (sqrt(2)*r/w1).^abs(l) ...
      .* laguerre_matrix(p, abs(l), 2*r.^2/w1^2) ...
      .* exp(-r.^2/w1^2);

% discrete normalisation
u = u / sqrt(sum(sum(conj(u).*u)));

% propagation phase: spiral + wavefront curvature + Gouy
phase = l*phi + (k/(2*R1))*r.^2 - (abs(l)+2*p+1)*atan(z1/zr);
u = u .* exp(1i*phase);
end
