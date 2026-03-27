function tl = pressure_to_tl(pressure)
%PRESSURE_TO_TL Convert complex pressure to transmission loss in dB.
tl = -20.0 * log10(abs(pressure) + eps);
end
