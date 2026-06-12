function I = add_gaussian_noise(I, snr_db)
noise   = abs(randn(size(I)));
var_n   = mean(noise(:).^2) - mean(noise(:))^2;
noise   = noise / sqrt(var_n);
var_s   = mean(I(:).^2) - mean(I(:))^2;
I       = I + (var_s / 10^(snr_db/10)) * noise;
I       = I / sum(sum(I));
end