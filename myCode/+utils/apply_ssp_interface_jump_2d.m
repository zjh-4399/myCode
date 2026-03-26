function yOut = apply_ssp_interface_jump_2d(yIn, hit, config)
%APPLY_SSP_INTERFACE_JUMP_2D Bellhop-style weak-interface jump for P.

yOut = yIn;

if ~local_jump_enabled(config)
    return;
end
if ~strcmpi(config.env.profile.type, 'tabular')
    return;
end
if ~strcmpi(local_get_tabular_mode(config), 'piecewise_linear')
    return;
end
if ~isfield(hit, 'interfaceIndex') || isempty(hit.interfaceIndex) || ~isfinite(hit.interfaceIndex)
    return;
end

zw = config.env.profile.zw(:);
cw = config.env.profile.cw(:);
j = hit.interfaceIndex;
if j <= 1 || j >= numel(zw)
    return;
end

r  = yIn(1);
z  = hit.zHit;
pr = yIn(3);
pz = yIn(4);
P  = complex(yIn(5), yIn(6));
Q  = complex(yIn(7), yIn(8));

if abs(pz) < 1.0e-12
    return;
end

dzUp = zw(j)   - zw(j - 1);
dzDn = zw(j+1) - zw(j);
if dzUp <= 0 || dzDn <= 0
    error('Invalid tabular depth grid around interface index %d.', j);
end

czAbove = (cw(j)   - cw(j - 1)) / dzUp;
czBelow = (cw(j+1) - cw(j    )) / dzDn;
[cIf, ~, ~, ~, ~, ~] = sound_speed_derivatives(config, r, z);

if pz > 0
    deltaCz = czBelow - czAbove;
else
    deltaCz = czAbove - czBelow;
end

RN = -(pr^2 / pz) * (deltaCz / cIf);
PNew = P + Q * RN;

yOut(5) = real(PNew);
yOut(6) = imag(PNew);
end

function tf = local_jump_enabled(config)
tf = true;
if isfield(config, 'numerics') && isfield(config.numerics, 'ssp_interface_jump') ...
        && ~isempty(config.numerics.ssp_interface_jump)
    tf = logical(config.numerics.ssp_interface_jump);
end
end

function mode = local_get_tabular_mode(config)
mode = 'piecewise_linear';
if isfield(config, 'env') && isfield(config.env, 'profile') ...
        && isfield(config.env.profile, 'tabular_mode') ...
        && ~isempty(config.env.profile.tabular_mode)
    mode = lower(strtrim(char(string(config.env.profile.tabular_mode))));
end
end
