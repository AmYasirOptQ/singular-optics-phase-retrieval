function [U1, meta] = random_lg_superposition(X, Y, dx, lambda, sup, Lin, Pin)
% RANDOM_LG_SUPERPOSITION  Random normalised superposition of LG modes.
%
%   [U1, meta] = random_lg_superposition(X, Y, dx, lambda, sup, Lin, Pin)
%
%   Inputs
%     X, Y   : N x N coordinate matrices (from make_grid)
%     dx     : grid spacing (m)
%     lambda : wavelength (m)
%     sup    : number of LG modes in the superposition
%     Lin    : 1 x M vector of available azimuthal indices l
%     Pin    : 1 x M vector of corresponding radial indices p
%
%   Outputs
%     U1   : N x N normalised complex superposition field
%     meta : struct with fields
%              .l_indices    – chosen azimuthal indices  (1 x sup)
%              .p_indices    – chosen radial indices     (1 x sup)
%              .coefficients – complex superposition coefficients (1 x sup)
%              .shifts       – transverse offsets in pixels (1 x 2*sup)
%              .w0           – beam waist (m)

N    = size(X, 1);
num0 = numel(Lin);

%% Random beam waist (Rayleigh range between 1e-4 and 3e-4 m)
w0     = 1e-4 + 2e-4*rand(1);
zR     = pi*w0^2/lambda;
z_beam = 0.001;
w1     = w0*sqrt(1+(z_beam/zR)^2);
R1     = z_beam*(1+(zR/z_beam)^2);

%% Random transverse offsets (up to 1% of N pixels per mode)
shift   = floor(0.01*N);
rand_no = floor(1 + (2*shift-1)*rand(1, 2*sup));

%% Randomly pick sup distinct LG modes
available = 1:num0;
Lind = zeros(1, sup);
Pind = zeros(1, sup);
for n = 1:sup
    idx      = floor(1 + (numel(available)-1)*rand(1));
    pick     = available(idx);
    Lind(n)  = Lin(pick);
    Pind(n)  = Pin(pick);
    available(available==pick) = [];
end

%% Generate each LG mode on its offset grid
LaF = zeros(N, N, sup);
for n = 1:sup
    Xs = X - rand_no(n)*dx;
    Ys = Y - rand_no(n+sup)*dx;
    rs   = sqrt(Xs.^2 + Ys.^2);
    phis = atan2(Ys, Xs);
    LaF(:,:,n) = LG_beam(rs, phis, Lind(n), Pind(n), w1, z_beam, zR, R1, lambda);
    LaF(:,:,n) = LaF(:,:,n) / sqrt(sum(sum(conj(LaF(:,:,n)).*LaF(:,:,n))));
end

%% Random complex superposition coefficients (unit sphere construction)
rc  = rand(1, sup-1)*pi;
rd  = exp(1i*rand(1, sup-1)*2*pi);
rd1 = [1, rd];

coeff = ones(sup-1, sup);
for va1 = 1:sup-1
    for va2 = 1:sup
        if va1==va2, coeff(va1,va2) = cos(rc(va1)); end
        if va1< va2, coeff(va1,va2) = sin(rc(va1)); end
    end
end

vect = ones(1, sup);
for va2 = 1:sup
    for va1 = 1:sup-1
        vect(va2) = vect(va2)*coeff(va1,va2);
    end
end
vect = vect.*rd1;

%% Form superposition
U1 = zeros(N);
for n = 1:sup
    U1 = U1 + LaF(:,:,n)*vect(n);
end
U1 = U1 / sqrt(sum(sum(conj(U1).*U1)));

%% Pack metadata
meta.l_indices    = Lind;
meta.p_indices    = Pind;
meta.coefficients = vect;
meta.shifts       = rand_no;
meta.w0           = w0;

end
