function L = beam_radius(config, P, Q)
%BEAM_RADIUS Return paraxial beam radius from P/Q.
pq = P ./ Q;
L = nan(size(pq));
mask = imag(pq) < 0;
L(mask) = sqrt(-2.0 ./ (config.omega .* imag(pq(mask))));
end
