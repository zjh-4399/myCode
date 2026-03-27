function [c, cr, cz, crr, crz, czz] = sound_speed_derivatives(config, r, z)
%SOUND_SPEED_DERIVATIVES
% Returns c and its derivatives for the current 2D solver.
%
% Supported profile types:
%   - downward
%   - munk
%   - tabular
%
% For tabular profiles, two modes are supported:
%   config.env.profile.tabular_mode = 'piecewise_linear'   % Bellhop-like weak interfaces
%   config.env.profile.tabular_mode = 'pp_spline'          % keep old ppval path
%
% In piecewise_linear mode:
%   c   : linear in each depth segment
%   cz  : piecewise constant
%   czz : zero inside each segment
% and the discontinuity of cz at internal interfaces is handled separately
% by apply_ssp_interface_jump_2d().

% 优化：用 persistent 缓存 profileType/tabMode/precomputed 标志，
% 消除每次调用重复执行的 lower()、isfield() 等固定开销（每步被调用 3-5 次）。
% 同一 MATLAB 会话中切换 config 时，可调用 clear sound_speed_derivatives 重置缓存。
persistent cacheProfileType cacheTabMode cacheUsePrecomputed
if isempty(cacheProfileType)
    cacheProfileType    = lower(config.env.profile.type);
    cacheTabMode        = local_get_tabular_mode(config);
    cacheUsePrecomputed = strcmp(cacheTabMode, 'piecewise_linear') && ...
        isfield(config.env.profile.params, 'pwl_zw')    && ...
        isfield(config.env.profile.params, 'pwl_cw')    && ...
        isfield(config.env.profile.params, 'pwl_slope') && ...
        isfield(config.env.profile.params, 'pwl_nseg');
end
profileType = cacheProfileType;
params = config.env.profile.params;

% 位置缓存：消除 local_predict_event 与积分器首步在相同 (r,z) 处的重复计算
persistent posR posZ posC posCr posCz posCrr posCrz posCzz
if isscalar(r) && isscalar(z) && ~isempty(posR) && r == posR && z == posZ
    c = posC; cr = posCr; cz = posCz;
    crr = posCrr; crz = posCrz; czz = posCzz;
    return;
end

switch profileType
    case 'downward'
        c0 = params.c0;
        gamma = params.gamma;
        arg = 1.0 - 2.0 * gamma .* z ./ c0;
        if any(arg <= 0)
            error('Downward profile invalid at z = %g.', z(1));
        end
        c = c0 ./ sqrt(arg);
        cz = gamma .* arg.^(-3/2);
        czz = (3.0 * gamma.^2 ./ c0) .* arg.^(-5/2);

    case 'munk'
        c0 = params.c0;
        eps_munk = params.epsilon;
        zc = params.zc;
        z_flat = params.z_flat;

        x = 2.0 * (z - zc) ./ zc;
        c = zeros(size(z));
        cz = zeros(size(z));
        czz = zeros(size(z));

        mask = z <= z_flat;
        c(mask) = c0 .* (1.0 + eps_munk .* (x(mask) - 1.0 + exp(-x(mask))));
        cz(mask) = c0 .* eps_munk .* (1.0 - exp(-x(mask))) .* (2.0 / zc);
        czz(mask) = c0 .* eps_munk .* exp(-x(mask)) .* (2.0 / zc)^2;

        if any(~mask)
            x_flat = 2.0 * (z_flat - zc) ./ zc;
            c_flat = c0 .* (1.0 + eps_munk .* (x_flat - 1.0 + exp(-x_flat)));
            c(~mask) = c_flat;
            cz(~mask) = 0.0;
            czz(~mask) = 0.0;
        end

    case 'tabular'
        % 优化：tabMode 和 precomputed 标志已在首次调用时缓存，无需重复检查
        if cacheUsePrecomputed
            [c, cz, czz] = local_tabular_piecewise_linear_precomputed( ...
                params.pwl_zw, params.pwl_cw, params.pwl_slope, params.pwl_nseg, z);
        elseif strcmp(cacheTabMode, 'piecewise_linear')
            [c, cz, czz] = local_tabular_piecewise_linear( ...
                config.env.profile.zw, config.env.profile.cw, z);
        elseif any(strcmp(cacheTabMode, {'pp_spline', 'spline_pp', 'pp'}))
            zQuery = min(max(z, min(config.env.profile.zw)), max(config.env.profile.zw));
            c   = ppval(params.ppc,   zQuery);
            cz  = ppval(params.ppcz,  zQuery);
            czz = ppval(params.ppczz, zQuery);
        else
            error('Unsupported tabular_mode: %s', cacheTabMode);
        end

    otherwise
        error('Unsupported env.profile.type: %s', profileType);
end

cr  = zeros(size(r));
crr = zeros(size(r));
crz = zeros(size(r));

% 更新位置缓存（仅标量查询）
if isscalar(r) && isscalar(z)
    posR = r; posZ = z;
    posC = c; posCr = cr; posCz = cz;
    posCrr = crr; posCrz = crz; posCzz = czz;
end
end

% =========================================================================
% Local helpers
% =========================================================================

function mode = local_get_tabular_mode(config)
mode = 'piecewise_linear';

if isfield(config, 'env') && isfield(config.env, 'profile') ...
        && isfield(config.env.profile, 'tabular_mode') ...
        && ~isempty(config.env.profile.tabular_mode)
    mode = lower(strtrim(char(string(config.env.profile.tabular_mode))));
end
end

function [c, cz, czz] = local_tabular_piecewise_linear(zw, cw, z)
zw = zw(:);
cw = cw(:);

if numel(zw) ~= numel(cw)
    error('env.profile.zw and env.profile.cw must have the same length.');
end
if numel(zw) < 2
    error('tabular profile must contain at least 2 depth points.');
end
if any(diff(zw) <= 0)
    error('env.profile.zw must be strictly increasing.');
end

zSize = size(z);
zCol = z(:);

if isscalar(zCol)
    zq = min(max(zCol, zw(1)), zw(end));
    if zq >= zw(end)
        segId = numel(zw) - 1;
    else
        segId = find(zw <= zq, 1, 'last');
        segId = min(max(segId, 1), numel(zw) - 1);
    end

    z0 = zw(segId);
    z1 = zw(segId + 1);
    c0 = cw(segId);
    c1 = cw(segId + 1);
    slope = (c1 - c0) / (z1 - z0);

    c = c0 + slope * (zq - z0);
    cz = slope;
    czz = 0.0;
    c = reshape(c, zSize);
    cz = reshape(cz, zSize);
    czz = reshape(czz, zSize);
    return;
end

% Clamp to valid interval. For exactly the last node, force it into the last segment.
zMin = zw(1);
zMax = zw(end);
zCol = min(max(zCol, zMin), zMax);

nSeg = numel(zw) - 1;

% discretize: bins are [zw(i), zw(i+1)) except the last edge
segId = discretize(zCol, zw);
maskLast = (zCol >= zMax);
segId(maskLast) = nSeg;

% For z == zw(1), discretize gives segId = 1 already
if any(isnan(segId))
    % very rare edge case due to floating-point equality at left edge
    segId(isnan(segId) & zCol <= zMin) = 1;
    segId(isnan(segId) & zCol >= zMax) = nSeg;
end

if any(isnan(segId))
    error('Failed to bracket some query depths in local_tabular_piecewise_linear.');
end

z0 = zw(segId);
z1 = zw(segId + 1);
c0 = cw(segId);
c1 = cw(segId + 1);

slope = (c1 - c0) ./ (z1 - z0);

cCol = c0 + slope .* (zCol - z0);
czCol = slope;
czzCol = zeros(size(cCol));

c   = reshape(cCol,  zSize);
cz  = reshape(czCol, zSize);
czz = reshape(czzCol, zSize);
end

function [c, cz, czz] = local_tabular_piecewise_linear_precomputed(zw, cw, slopeSeg, nSeg, z)
zSize = size(z);
zCol = z(:);

zMin = zw(1);
zMax = zw(end);
zCol = min(max(zCol, zMin), zMax);

if isscalar(zCol)
    if zCol >= zMax
        segId = nSeg;
    else
        segId = local_cached_segment_id(zw, nSeg, zCol);
    end

    cVal = cw(segId) + slopeSeg(segId) * (zCol - zw(segId));
    c = reshape(cVal, zSize);
    cz = reshape(slopeSeg(segId), zSize);
    czz = reshape(0.0, zSize);
    return;
end

segId = discretize(zCol, zw);
segId(zCol >= zMax) = nSeg;
segId(isnan(segId) & zCol <= zMin) = 1;
segId(isnan(segId) & zCol >= zMax) = nSeg;

if any(isnan(segId))
    segId(isnan(segId)) = 1;
end

cCol = cw(segId) + slopeSeg(segId) .* (zCol - zw(segId));
czCol = slopeSeg(segId);
czzCol = zeros(size(cCol));

c = reshape(cCol, zSize);
cz = reshape(czCol, zSize);
czz = reshape(czzCol, zSize);
end

function segId = local_cached_segment_id(zw, nSeg, zq)
% Bellhop-like layer chasing: reuse the previous segment and move up/down.
persistent lastSeg
if isempty(lastSeg) || ~isfinite(lastSeg)
    lastSeg = 1;
end
if lastSeg < 1 || lastSeg > nSeg
    lastSeg = 1;
end

while (lastSeg < nSeg) && (zq >= zw(lastSeg + 1))
    lastSeg = lastSeg + 1;
end
while (lastSeg > 1) && (zq < zw(lastSeg))
    lastSeg = lastSeg - 1;
end

segId = lastSeg;
end
