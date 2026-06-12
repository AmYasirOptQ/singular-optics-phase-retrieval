%% main.m
%
% Three-plane phase retrieval using Laguerre-Gaussian beams,
% a partial Fourier transform, and Fresnel propagation.
%
% Plane convention (as in the paper):
%
%   LG superposition  (plane 1)
%      |
%      |  Partial Fourier Transform
%      v
%   Rotated HG-like field  (plane 2)
%      |
%      |  Fresnel propagation
%      v
%   Propagated field  (plane 3)
%      |
%      +--> Intensities I1, I2, I3 (+ noise + aperture)
%      |
%      +--> Gerchberg-Saxton reconstruction
%      |
%      +--> Twist parameter comparison

clear; clc;

%% ------------------------------------------------------------------------
% Simulation parameters
%% ------------------------------------------------------------------------

N      = 512;
iter   = 1;       % number of independent random samples
sup    = 2;       % number of LG modes in superposition
num    = 5;       % highest j-index: j in {0, 1/2, ..., num/2}
count  = 100;     % GS iterations per sample

lambda      = 632.8e-9;
focal       = 0.1*sqrt(2);
L           = sqrt(N*lambda*focal);
z           = focal;
snr_db      = 10;
mask_radius = floor(0.2*N);

%% ------------------------------------------------------------------------
% Computational grid
%% ------------------------------------------------------------------------

[X, Y, r, phi, dx] = make_grid(N, L);

%% ------------------------------------------------------------------------
% LG index table:  (j, m) -> (l, p)
%% ------------------------------------------------------------------------

[Lin, Pin] = lg_mode_table(num);

%% ------------------------------------------------------------------------
% Pre-allocate output arrays
%% ------------------------------------------------------------------------

chin1=zeros(1,iter);  chin2=zeros(1,iter);  chin3=zeros(1,iter);
chout1=zeros(1,iter); chout2=zeros(1,iter); chout3=zeros(1,iter);
error1=zeros(1,iter); error2=zeros(1,iter); error3=zeros(1,iter);
inten1=zeros(1,iter); inten2=zeros(1,iter); inten3=zeros(1,iter);
ampl1=zeros(1,iter);  ampl2=zeros(1,iter);  ampl3=zeros(1,iter);
col1=zeros(1,iter);
eflag = 0;

%% ------------------------------------------------------------------------
% Main loop over independent samples
%% ------------------------------------------------------------------------

for it = 1:iter

    disp('iteration number'); disp(it);

    %% --------------------------------------------------------------------
    % Plane 1: random LG superposition
    %% --------------------------------------------------------------------

    [U1, meta] = random_lg_superposition(X, Y, dx, lambda, sup, Lin, Pin);

    %% --------------------------------------------------------------------
    % Plane 2: partial Fourier transform
    % Plane 3: Fresnel propagation
    %% --------------------------------------------------------------------

    U2 = PFTransform(U1, N);
    U3 = fresnel_propagation(U2, z, lambda, L);
    U3 = U3 / sqrt(sum(sum(conj(U3).*U3)));

    %% --------------------------------------------------------------------
    % Input twist parameters
    %% --------------------------------------------------------------------

    tau01 = twist_parameter(U1, X, Y, dx, lambda);
    tau02 = twist_parameter(U2, X, Y, dx, lambda);
    tau03 = twist_parameter(U3, X, Y, dx, lambda);

    %% --------------------------------------------------------------------
    % Measured intensities
    %% --------------------------------------------------------------------

    I1 = abs(U1).^2;
    I2 = abs(U2).^2;
    I3 = abs(U3).^2;

    %% --------------------------------------------------------------------
    % Add Gaussian noise
    %% --------------------------------------------------------------------

    I1 = add_gaussian_noise(I1, snr_db);
    I2 = add_gaussian_noise(I2, snr_db);
    I3 = add_gaussian_noise(I3, snr_db);

    %% --------------------------------------------------------------------
    % Apply circular aperture
    %% --------------------------------------------------------------------

    mask = circ_aper(N, mask_radius);

    I1 = normalize_intensity(I1 .* mask);
    I2 = normalize_intensity(I2 .* mask);
    I3 = normalize_intensity(I3 .* mask);

    U1t = U1.*mask;  U1t = U1t/sqrt(sum(sum(conj(U1t).*U1t)));
    U2t = U2.*mask;  U2t = U2t/sqrt(sum(sum(conj(U2t).*U2t)));
    U3t = U3.*mask;  U3t = U3t/sqrt(sum(sum(conj(U3t).*U3t)));

    %% --------------------------------------------------------------------
    % Define operators
    %% --------------------------------------------------------------------

    PF        = @(U) PFTransform(U, N);
    PF_inv    = @(U) IPFTransform(U, N);
    propagate = @(U) fresnel_propagation(U,  z, lambda, L);
    backprop  = @(U) fresnel_propagation(U, -z, lambda, L);

    %% --------------------------------------------------------------------
    % Gerchberg-Saxton reconstruction
    %% --------------------------------------------------------------------

    options.display = true;

    [U_rec, history] = gerchberg_saxton_3plane( ...
        I1, I2, I3, ...
        PF, PF_inv, ...
        propagate, backprop, ...
        count, options);

    %% --------------------------------------------------------------------
    % Twist parameters of reconstructed fields
    %% --------------------------------------------------------------------

    tau1 = twist_parameter(U_rec, X, Y, dx, lambda);

    U2_rec = PFTransform(U_rec, N);
    U3_rec = fresnel_propagation(U2_rec, z, lambda, L);
    tau2   = twist_parameter(U2_rec, X, Y, dx, lambda);
    tau3   = twist_parameter(U3_rec, X, Y, dx, lambda);

    %% --------------------------------------------------------------------
    % Store results
    %% --------------------------------------------------------------------

    chin1(it)=tau01;  chin2(it)=tau02;  chin3(it)=tau03;
    chout1(it)=tau1;  chout2(it)=tau2;  chout3(it)=tau3;

    error1(it) = abs(tau01-tau1)*100/abs(tau01);
    error2(it) = abs(tau02-tau2)*100/abs(tau02);
    error3(it) = abs(tau03-tau3)*100/abs(tau03);

    inten1(it) = history.corr_plane1(count);
    inten2(it) = history.corr_plane2(count);
    inten3(it) = history.corr_plane3(count);

    ampl1(it)  = history.ampl_plane1(count);
    ampl2(it)  = history.ampl_plane2(count);
    ampl3(it)  = history.ampl_plane3(count);

    col1(it) = abs(chin1(it)-chout1(it));

    if (col1(it)/abs(chin1(it))*100) > 5
        eflag=eflag+1;
        eMeta(eflag) = meta; %#ok<SAGROW>
    end

    %% --------------------------------------------------------------------
    % Report
    %% --------------------------------------------------------------------

    fprintf('\n--------------------------------------\n');
    fprintf('Sample %d / %d\n', it, iter);
    fprintf('l indices : %s\n', num2str(meta.l_indices));
    fprintf('p indices : %s\n', num2str(meta.p_indices));
    fprintf('Twist parameter input  (plane 1): %10.6f\n', real(tau01));
    fprintf('Twist parameter recon  (plane 1): %10.6f\n', real(tau1));
    fprintf('Relative error                  : %10.4f %%\n', error1(it));
    fprintf('--------------------------------------\n');

end

%% ------------------------------------------------------------------------
% Summary
%% ------------------------------------------------------------------------

fprintf('\nMean intensity corr (plane 1): %.4f\n', mean(inten1));
fprintf('Mean intensity corr (plane 2): %.4f\n', mean(inten2));
fprintf('Mean intensity corr (plane 3): %.4f\n', mean(inten3));
fprintf('Mean charge error%% (plane 1):  %.2f\n',  mean(error1));
fprintf('Non-converged samples: %d / %d\n', eflag, iter);

%% ------------------------------------------------------------------------
% Visualization
%% ------------------------------------------------------------------------

figure;
subplot(2,2,1); imagesc(abs(U1).^2);    axis image; colorbar; title('Input intensity (plane 1)');
subplot(2,2,2); imagesc(angle(U1));     axis image; colorbar; title('Input phase (plane 1)');
subplot(2,2,3); imagesc(abs(U_rec).^2); axis image; colorbar; title('Reconstructed intensity (plane 1)');
subplot(2,2,4); imagesc(angle(U_rec));  axis image; colorbar; title('Reconstructed phase (plane 1)');

figure;
plot(1:count, history.corr_plane1, 'LineWidth', 1.5); hold on;
plot(1:count, history.corr_plane2, 'LineWidth', 1.5);
plot(1:count, history.corr_plane3, 'LineWidth', 1.5);
xlabel('Iteration'); ylabel('Intensity correlation');
legend('Plane 1 (LG)', 'Plane 2 (PF)', 'Plane 3 (propagated)');
title('GS convergence'); grid on;

