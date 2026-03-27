function yOut = apply_ssp_interface_jump_2d(yIn, hit, config)
%APPLY_SSP_INTERFACE_JUMP_2D Bellhop-style weak-interface jump for P.

yOut = yIn;

% 优化：缓存 jump 启用标志及 zw/cw 数组，消除每次界面穿越的 isfield/strcmpi 检查与数组拷贝
% 切换 config 时可调用 clear utils.apply_ssp_interface_jump_2d 重置缓存
persistent jumpApplicable cachedZw cachedCw cachedNZw
if isempty(jumpApplicable)
    jumpApplicable = local_jump_enabled(config) && ...
        strcmpi(config.env.profile.type, 'tabular') && ...
        strcmpi(local_get_tabular_mode(config), 'piecewise_linear');
    if jumpApplicable
        cachedZw  = config.env.profile.zw(:);
        cachedCw  = config.env.profile.cw(:);
        cachedNZw = numel(cachedZw);
    else
        cachedZw  = [];
        cachedCw  = [];
        cachedNZw = 0;
    end
end
if ~jumpApplicable
    return;
end
if ~isfield(hit, 'interfaceIndex') || isempty(hit.interfaceIndex) || ~isfinite(hit.interfaceIndex)
    return;
end

zw = cachedZw;
cw = cachedCw;
j = hit.interfaceIndex;
if j <= 1 || j >= cachedNZw
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
