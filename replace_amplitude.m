function u = replace_amplitude(u, I)
% Keep phase, replace amplitude with sqrt(I).
u = exp(1i * atan2(imag(u), real(u))) .* sqrt(I);
end