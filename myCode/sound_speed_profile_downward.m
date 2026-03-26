function c = sound_speed_profile_downward(z, params)
%SOUND_SPEED_PROFILE_DOWNWARD Downward refracting profile.
c0 = params.c0;
gamma = params.gamma;
arg = 1.0 - 2.0 * gamma .* z ./ c0;
if any(arg <= 0)
    error('Downward profile invalid for the supplied z range.');
end
c = c0 ./ sqrt(arg);
end
