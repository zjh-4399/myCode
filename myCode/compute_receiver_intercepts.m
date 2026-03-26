function intercepts = compute_receiver_intercepts(beam, receiver_depth)
%COMPUTE_RECEIVER_INTERCEPTS Compute where beam normals hit z = receiver_depth.
Nz = beam.Nz;
valid = isfinite(Nz) & abs(Nz) > 1.0e-10;
r_int = nan(size(beam.r));
n_int = nan(size(beam.r));
r_int(valid) = beam.r(valid) + (receiver_depth - beam.z(valid)) .* beam.Nr(valid) ./ Nz(valid);
n_int(valid) = (receiver_depth - beam.z(valid)) ./ Nz(valid);

intercepts = struct();
intercepts.valid = valid;
intercepts.r_int = r_int;
intercepts.n_int = n_int;
end
