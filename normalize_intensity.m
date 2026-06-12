function I = normalize_intensity(I)
s = sum(sum(I));
if s > 0, I = I / s; end
end