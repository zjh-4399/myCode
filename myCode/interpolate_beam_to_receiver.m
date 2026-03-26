function [uLine, debugInfo] = interpolate_beam_to_receiver(beam, receiver_depth, receiver_ranges, config)
%INTERPOLATE_BEAM_TO_RECEIVER Interpolate one beam to a receiver line.
%
% 加速版思路：
% 1) 对均匀 receiver_ranges，用索引映射代替每段里的 find(...)
% 2) beam/config 元数据只解析一次
% 3) 段内插值先向量化，再逐点做 branch tracking / contribution
%
% 说明：
% - 保持与原版相同的物理公式与 continuous_sqrt 分支连续跟踪逻辑
% - 非均匀 receiver_ranges 自动回退到通用查找
% - 可直接替换原 interpolate_beam_to_receiver.m

nR = numel(receiver_ranges);
uLine = complex(zeros(1, nR));

debugInfo = struct();
debugInfo.receiver_depth = receiver_depth;
debugInfo.segments_used  = 0;
debugInfo.points_used    = 0;

if nR == 0 || numel(beam.r) < 2
    return;
end

% ----------------------------
% 预处理：接收线截距、beam 元数据、receiver 网格信息
% ----------------------------
intercepts = compute_receiver_intercepts_local(beam, receiver_depth);
beamMeta   = local_unpack_beam_meta(config);
gridInfo   = local_receiver_grid_info(receiver_ranges);

lastRoot = [];
segmentsUsed = 0;
pointsUsed   = 0;

nSeg = numel(beam.r) - 1;

for j = 1:nSeg
    if ~(intercepts.valid(j) && intercepts.valid(j+1))
        continue;
    end

    ra = intercepts.r_int(j);
    rb = intercepts.r_int(j+1);

    if ~isfinite(ra) || ~isfinite(rb)
        continue;
    end

    drab = rb - ra;
    if abs(drab) < 1.0e-12
        continue;
    end

    rMin = min(ra, rb);
    rMax = max(ra, rb);

    % ----------------------------
    % 关键加速：均匀网格直接算索引区间
    % ----------------------------
    [i1, i2] = local_receiver_index_span(receiver_ranges, gridInfo, rMin, rMax);
    if i2 < i1
        continue;
    end

    idx = i1:i2;
    rr  = receiver_ranges(idx);

    % 段内线性插值（向量化）
    w = (rr - ra) ./ drab;

    G     = (1.0 - w) .* beam.gain(j) + w .* beam.gain(j+1);
    n     = (1.0 - w) .* intercepts.n_int(j) + w .* intercepts.n_int(j+1);
    P     = (1.0 - w) .* beam.P(j) + w .* beam.P(j+1);
    Q     = (1.0 - w) .* beam.Q(j) + w .* beam.Q(j+1);
    tau   = (1.0 - w) .* beam.tau(j) + w .* beam.tau(j+1);
    c     = (1.0 - w) .* beam.c(j) + w .* beam.c(j+1);
    rBeam = (1.0 - w) .* beam.r(j) + w .* beam.r(j+1);

    segmentsUsed = segmentsUsed + 1;

    % ----------------------------
    % 逐点贡献（保留 continuous_sqrt 的分支连续跟踪）
    % ----------------------------
    for k = 1:numel(idx)
        [contrib, lastRoot, used] = local_beam_contribution_fast( ...
            beamMeta, config, beam, ...
            n(k), P(k), Q(k), tau(k), c(k), rBeam(k), lastRoot, G(k));

        if used
            uLine(idx(k)) = uLine(idx(k)) + contrib;
            pointsUsed = pointsUsed + 1;
        end
    end
end

debugInfo.segments_used = segmentsUsed;
debugInfo.points_used   = pointsUsed;

end


% =========================================================================
%                           Local helpers
% =========================================================================

function intercepts = compute_receiver_intercepts_local(beam, receiver_depth)
% 计算 beam 法线与 z = receiver_depth 的截距
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


function meta = local_unpack_beam_meta(config)
% 只解析一次 beam 元数据，避免每个接收点反复 isfield 判断

meta = struct();

% family
if isfield(config, 'beam') && isfield(config.beam, 'model') && isfield(config.beam.model, 'family')
    meta.family = lower(char(string(config.beam.model.family)));
elseif isfield(config, 'beam') && isfield(config.beam, 'family')
    meta.family = lower(char(string(config.beam.family)));
else
    error('Cannot find beam family in config.beam.family or config.beam.model.family.');
end

% 统一读字段：优先 model，其次 beam
meta.autoWindow = local_read_beam_field(config, 'auto_window', false);
meta.windowRadii = local_read_beam_field(config, 'window_radii', 5.0);
meta.widthScale  = local_read_beam_field(config, 'width_scale', 1.0);
meta.widthFloor  = local_read_beam_field(config, 'width_floor', 0.0);
meta.useStent    = local_read_beam_field(config, 'use_stent', false);

if isfield(config.beam, 'delta_alpha_rad') && ~isempty(config.beam.delta_alpha_rad)
    meta.deltaAlpha = config.beam.delta_alpha_rad;
else
    error('config.beam.delta_alpha_rad is missing. Make sure validate_config_2d has run.');
end
end


function value = local_read_beam_field(config, fieldName, defaultValue)
if isfield(config, 'beam') && isfield(config.beam, 'model') && isfield(config.beam.model, fieldName) ...
        && ~isempty(config.beam.model.(fieldName))
    value = config.beam.model.(fieldName);
elseif isfield(config, 'beam') && isfield(config.beam, fieldName) ...
        && ~isempty(config.beam.(fieldName))
    value = config.beam.(fieldName);
else
    value = defaultValue;
end
end


function gridInfo = local_receiver_grid_info(receiver_ranges)
% 判断是否为严格递增、近似均匀网格
gridInfo = struct();
gridInfo.isUniform = false;
gridInfo.isIncreasing = false;
gridInfo.r0 = [];
gridInfo.dr = [];

rr = receiver_ranges(:);
if numel(rr) <= 1
    gridInfo.isUniform = true;
    gridInfo.isIncreasing = true;
    if ~isempty(rr)
        gridInfo.r0 = rr(1);
    else
        gridInfo.r0 = 0.0;
    end
    gridInfo.dr = 1.0;
    return;
end

d = diff(rr);
gridInfo.isIncreasing = all(d > 0);

if ~gridInfo.isIncreasing
    return;
end

d0 = d(1);
tol = 1.0e-10 * max(1.0, abs(d0));
gridInfo.isUniform = all(abs(d - d0) <= tol);
gridInfo.r0 = rr(1);
gridInfo.dr = d0;
end


function [i1, i2] = local_receiver_index_span(receiver_ranges, gridInfo, rMin, rMax)
% 返回 receiver_ranges 中落在 [rMin, rMax] 的索引区间 [i1, i2]
nR = numel(receiver_ranges);

if nR == 0
    i1 = 1;
    i2 = 0;
    return;
end

if gridInfo.isUniform && gridInfo.isIncreasing
    r0 = gridInfo.r0;
    dr = gridInfo.dr;

    x1 = (rMin - r0) / dr + 1.0;
    x2 = (rMax - r0) / dr + 1.0;

    % 加一点 tiny 容差，避免端点因浮点误差被漏掉
    i1 = max(1, ceil (x1 - 1.0e-12));
    i2 = min(nR, floor(x2 + 1.0e-12));

    if i2 < i1
        i1 = 1;
        i2 = 0;
    end
else
    % 非均匀网格：回退到通用查找
    i1 = find(receiver_ranges >= rMin, 1, 'first');
    i2 = find(receiver_ranges <= rMax, 1, 'last');

    if isempty(i1) || isempty(i2) || i2 < i1
        i1 = 1;
        i2 = 0;
    end
end
end


function [contrib, lastRoot, used] = local_beam_contribution_fast( ...
    meta, config, beam, n, P, Q, tau, c, rBeam, lastRoot, G)

contrib = 0;
used = false;

if ~isfinite(rBeam) || rBeam <= 0 || ~isfinite(c) || ~isfinite(P) || ~isfinite(Q)
    return;
end

switch meta.family
    case 'paraxial'
        pq = P / Q;

        if ~(imag(pq) < 0)
            return;
        end

        L = sqrt(-2.0 / (config.omega * imag(pq)));
        if ~isfinite(L)
            return;
        end

        if meta.autoWindow
            if abs(n) > meta.windowRadii * L
                return;
            end
        end

        [rootVal, lastRoot] = utils.continuous_sqrt(c / (rBeam * Q), lastRoot);
        phase   = exp(-1i * config.omega * (tau + 0.5 * pq * n^2));
        contrib = G * meta.deltaAlpha * beam.Aconst * rootVal * phase;
        used    = true;

    case 'geometric_hat'
        [qEff, W] = local_geometric_width_fast(meta, real(Q), c, false);
        if ~isfinite(W) || W <= 0 || abs(n) > W
            return;
        end

        shape = (W - abs(n)) / W;
        [rootVal, lastRoot] = utils.continuous_sqrt(c / (rBeam * complex(qEff, 0.0)), lastRoot);
        amp     = sqrt(max(cos(beam.alpha_rad), 1.0e-12)) * rootVal;
        phase   = exp(-1i * config.omega * tau);
        contrib = G * beam.Aconst * amp * shape * phase;
        used    = true;

    case 'geometric_gaussian'
        [qEff, W] = local_geometric_width_fast(meta, real(Q), c, true);
        if ~isfinite(W) || W <= 0
            return;
        end

        shape = exp(-0.5 * (n / W)^2);
        [rootVal, lastRoot] = utils.continuous_sqrt(c / (rBeam * complex(qEff, 0.0)), lastRoot);
        amp     = (1.0 / sqrt(2.0 * pi)) * sqrt(max(cos(beam.alpha_rad), 1.0e-12)) * rootVal;
        phase   = exp(-1i * config.omega * tau);
        contrib = G * beam.Aconst * amp * shape * phase;
        used    = true;

    otherwise
        error('Unsupported beam family: %s', meta.family);
end
end


function [qEff, W] = local_geometric_width_fast(meta, qReal, c0local, useStentRequested)
dalpha = abs(meta.deltaAlpha);
if dalpha <= 0
    dalpha = deg2rad(1.0);
end

qAbs  = abs(qReal);
Wgeom = meta.widthScale * qAbs * dalpha / max(abs(c0local), eps);

if useStentRequested && meta.useStent
    W = max(Wgeom, meta.widthFloor);
    qFloor = max(abs(c0local) * W / max(meta.widthScale * dalpha, eps), eps);

    if qReal < 0
        qEff = -max(qAbs, qFloor);
    else
        qEff =  max(qAbs, qFloor);
    end
else
    W = Wgeom;
    if qReal == 0
        qEff = eps;
    else
        qEff = qReal;
    end
end
end