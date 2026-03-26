function c = sound_speed_profile_munk(z, params)
%SOUND_SPEED_PROFILE_MUNK Munk profile with flattening below z_flat.
c0 = params.c0;
eps_munk = params.epsilon;
zc = params.zc;
z_flat = params.z_flat;

x = 2.0 * (z - zc) ./ zc;
c = zeros(size(z));

mask = z <= z_flat;
c(mask) = c0 .* (1.0 + eps_munk .* (x(mask) - 1.0 + exp(-x(mask))));
if any(~mask)
    x_flat = 2.0 * (z_flat - zc) ./ zc;
    c_flat = c0 .* (1.0 + eps_munk .* (x_flat - 1.0 + exp(-x_flat)));
    c(~mask) = c_flat;
end
end
