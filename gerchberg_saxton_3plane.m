function [U_rec, history] = gerchberg_saxton_3plane( ...
    I1, I2, I3, PF, PF_inv, propagate, backprop, num_iterations, options)
% GERCHBERG_SAXTON_THREE_PLANE  Three-plane GS phase retrieval.
%
%   [U_rec, history] = gerchberg_saxton_three_plane( ...
%       I1, I2, I3, PF, PF_inv, propagate, backprop, num_iterations, options)
%
%   Plane convention (as in the paper):
%     Plane 1 : LG field          U1          intensity I1
%     Plane 2 : partial FT        U2=PF(U1)   intensity I2
%     Plane 3 : Fresnel prop.     U3=prop(U2) intensity I3
%
%   Each GS iteration:
%     1. U2  = PF(U)            -> replace amplitude with sqrt(I2)
%     2. U3  = propagate(U2)    -> replace amplitude with sqrt(I3)
%     3. U2b = backprop(U3)     -> replace amplitude with sqrt(I2)
%     4. U1  = PF_inv(U2b)      -> replace amplitude with sqrt(I1)
%
%   Inputs
%     I1, I2, I3     : N x N measured (apertured, normalised) intensities
%     PF             : function handle — forward partial Fourier transform
%     PF_inv         : function handle — inverse partial Fourier transform
%     propagate      : function handle — forward Fresnel propagation
%     backprop       : function handle — backward Fresnel propagation
%     num_iterations : number of GS iterations
%     options        : struct; options.display = true shows live phase plot
%
%   Outputs
%     U_rec   : N x N reconstructed complex field at plane 1
%     history : struct with fields
%                 .corr_plane1  – intensity correlation at plane 1 (LG)
%                 .corr_plane2  – intensity correlation at plane 2 (PF)
%                 .corr_plane3  – intensity correlation at plane 3 (propagated)
%                 .ampl_plane1  – amplitude overlap at plane 1
%                 .ampl_plane2  – amplitude overlap at plane 2
%                 .ampl_plane3  – amplitude overlap at plane 3

if nargin < 9 || ~isfield(options, 'display')
    options.display = false;
end

history.corr_plane1 = zeros(1, num_iterations);
history.corr_plane2 = zeros(1, num_iterations);
history.corr_plane3 = zeros(1, num_iterations);
history.ampl_plane1 = zeros(1, num_iterations);
history.ampl_plane2 = zeros(1, num_iterations);
history.ampl_plane3 = zeros(1, num_iterations);

% initialise with measured amplitude at plane 1, zero phase
U = sqrt(I1);

for m = 1 : num_iterations

    % --- plane 1 -> plane 2 : partial Fourier transform ---
    U2 = PF(U);
    history.corr_plane2(m) = corr2(I2, abs(U2).^2);
    history.ampl_plane2(m) = abs(sum(sum(U2 .* conj(sqrt(I2)))));
    U2 = replace_amplitude(U2, I2);

    % --- plane 2 -> plane 3 : Fresnel propagation ---
    U3 = propagate(U2);
    history.corr_plane3(m) = corr2(I3, abs(U3).^2);
    history.ampl_plane3(m) = abs(sum(sum(U3 .* conj(sqrt(I3)))));
    U3 = replace_amplitude(U3, I3);

    % --- plane 3 -> plane 2 : backward Fresnel propagation ---
    U2b = backprop(U3);
    U2b = replace_amplitude(U2b, I2);

    % --- plane 2 -> plane 1 : inverse partial Fourier transform ---
    U1  = PF_inv(U2b);
    history.corr_plane1(m) = corr2(I1, abs(U1).^2);
    history.ampl_plane1(m) = abs(sum(sum(U1 .* conj(sqrt(I1)))));
    U1  = replace_amplitude(U1, I1);

    U = U1;

    if options.display
        imagesc(angle(U));
        title(sprintf('GS iteration %d', m));
        drawnow;
    end

end

U_rec = U;
end

