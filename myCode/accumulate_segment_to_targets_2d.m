function targets = accumulate_segment_to_targets_2d(targets, pointA, pointB, beamInfo, config)
%ACCUMULATE_SEGMENT_TO_TARGETS_2D Add one traced ray segment directly to all targets.
%
% This is the segment-level online replacement for the older workflow:
%   trace full beam -> interpolate_beam_to_receiver over all segments.
%
% Here each newly generated segment [pointA, pointB] contributes
% immediately to the pressure accumulator.

if nargin < 1 || isempty(targets)
    return;
end

if ~isstruct(targets)
    error('targets must be a struct array returned by init_segment_accumulator_2d.');
end

for it = 1:numel(targets)
    targets(it) = local_accumulate_one_target(targets(it), pointA, pointB, beamInfo, config);
end
end


function target = local_accumulate_one_target(target, pointA, pointB, beamInfo, config)
meta = target.beamMeta;
rr = target.ranges;
nR = numel(rr);
if nR == 0
    return;
end

validA = isfinite(pointA.Nz) && abs(pointA.Nz) > 1.0e-10;
validB = isfinite(pointB.Nz) && abs(pointB.Nz) > 1.0e-10;
if ~(validA && validB)
    return;
end

depths = target.depths;
raAll = pointA.r + (depths - pointA.z) * pointA.Nr / pointA.Nz;
rbAll = pointB.r + (depths - pointB.z) * pointB.Nr / pointB.Nz;
naAll = (depths - pointA.z) / pointA.Nz;
nbAll = (depths - pointB.z) / pointB.Nz;

maskFinite = isfinite(raAll) & isfinite(rbAll) & isfinite(naAll) & isfinite(nbAll);
if ~any(maskFinite)
    return;
end

drabAll = rbAll - raAll;
rMinAll = min(raAll, rbAll);
rMaxAll = max(raAll, rbAll);
% 优化：将 drab 零检查合并进 maskOverlap，消除循环内每深度的分支判断
maskOverlap = maskFinite & (abs(drabAll) >= 1.0e-12) & ...
    (rMaxAll >= rr(1)) & (rMinAll <= rr(end));
if ~any(maskOverlap)
    return;
end

depthIdx = find(maskOverlap);
nActive = numel(depthIdx);

% 优化：对均匀接收距离网格，用 ceil/floor 算术一次性算出所有活跃深度的索引区间，
% 消除循环内逐深度调用 local_receiver_index_span 的开销
if target.gridInfo.isUniform && target.gridInfo.isIncreasing
    r0_g = target.gridInfo.r0;
    dr_g = target.gridInfo.dr;
    x1v = (rMinAll(depthIdx) - r0_g) / dr_g + 1.0;
    x2v = (rMaxAll(depthIdx) - r0_g) / dr_g + 1.0;
    i1v = max(1, ceil(x1v - 1.0e-12));
    i2v = min(nR, floor(x2v + 1.0e-12));
% else
%     i1v = zeros(nActive, 1);
%     i2v = zeros(nActive, 1);
%     for jj = 1:nActive
%         iz = depthIdx(jj);
%         [i1v(jj), i2v(jj)] = local_receiver_index_span(rr, target.gridInfo, rMinAll(iz), rMaxAll(iz));
%     end
end
validSpan = (i2v >= i1v);
depthIdx  = depthIdx(validSpan);
i1v = i1v(validSpan);
i2v = i2v(validSpan);
if isempty(depthIdx)
    return;
end

% 优化：预计算段内常量差值，循环内改用 lerp 形式 A + w*(B-A)，
% 每行从 2 次 nR 向量乘法减为 1 次，同时将 12 次结构体字段访问移出循环
gainA      = pointA.gain;  gainDelta  = pointB.gain - pointA.gain;
PA         = pointA.P;     PDelta     = pointB.P    - pointA.P;
QA         = pointA.Q;     QDelta     = pointB.Q    - pointA.Q;
tauA       = pointA.tau;   tauDelta   = pointB.tau  - pointA.tau;
cA         = pointA.c;     cDelta     = pointB.c    - pointA.c;
rBeamA     = pointA.r;     rBeamDelta = pointB.r    - pointA.r;
% 供等比数列相位使用的常量（仅 hat/gaussian 分支需要）
omega       = config.omega;
isUniform_g = target.gridInfo.isUniform && target.gridInfo.isIncreasing;
dr_g        = target.gridInfo.dr;

% 优化：将 switch 提到循环外，每次 accumulate 调用只判断一次 familyCode
switch meta.familyCode
    case 1
        for ii = 1:numel(depthIdx)
            iz = depthIdx(ii);
            ra = raAll(iz);
            na = naAll(iz);
            nb = nbAll(iz);
            drab = drabAll(iz);
            idx = i1v(ii):i2v(ii);
            w = (rr(idx) - ra) ./ drab;
            G     = gainA  + w .* gainDelta;
            n     = na     + w .* (nb - na);
            P     = PA     + w .* PDelta;
            Q     = QA     + w .* QDelta;
            tau   = tauA   + w .* tauDelta;
            c     = cA     + w .* cDelta;
            rBeam = rBeamA + w .* rBeamDelta;
            lastRoot = target.lastRoots{iz};
            [linePressure, lastRoot] = local_accumulate_paraxial_line( ...
                target.pressure(iz, idx), meta, config, beamInfo, ...
                n, P, Q, tau, c, rBeam, lastRoot, G);
            target.pressure(iz, idx) = linePressure;
            target.lastRoots{iz} = lastRoot;
        end

    case 2
        for ii = 1:numel(depthIdx)
            iz = depthIdx(ii);
            ra = raAll(iz);
            na = naAll(iz);
            nb = nbAll(iz);
            drab = drabAll(iz);
            idx = i1v(ii):i2v(ii);
            w = (rr(idx) - ra) ./ drab;
            G     = gainA  + w .* gainDelta;
            n     = na     + w .* (nb - na);
            Q     = QA     + w .* QDelta;
            tau   = tauA   + w .* tauDelta;
            c     = cA     + w .* cDelta;
            rBeam = rBeamA + w .* rBeamDelta;
            % 优化：均匀网格用等比数列相位，将 nRange 次 exp 压缩为 2 次 exp + nRange 次复数乘法
            if isUniform_g
                nR_loc = numel(idx);
                phaseStart = exp(-1i * omega * tau(1));
                if nR_loc > 1
                    phaseStep = exp(-1i * omega * (dr_g / drab * tauDelta));
                    phase_g = cumprod([phaseStart, repmat(phaseStep, 1, nR_loc - 1)]);
                else
                    phase_g = phaseStart;
                end
            else
                phase_g = exp(-1i * omega * tau);
            end
            lastRoot = target.lastRoots{iz};
            [linePressure, lastRoot] = local_accumulate_geometric_hat_line( ...
                target.pressure(iz, idx), meta, beamInfo, ...
                n, Q, phase_g, c, rBeam, lastRoot, G);
            target.pressure(iz, idx) = linePressure;
            target.lastRoots{iz} = lastRoot;
        end

    case 3
        for ii = 1:numel(depthIdx)
            iz = depthIdx(ii);
            ra = raAll(iz);
            na = naAll(iz);
            nb = nbAll(iz);
            drab = drabAll(iz);
            idx = i1v(ii):i2v(ii);
            w = (rr(idx) - ra) ./ drab;
            G     = gainA  + w .* gainDelta;
            n     = na     + w .* (nb - na);
            Q     = QA     + w .* QDelta;
            tau   = tauA   + w .* tauDelta;
            c     = cA     + w .* cDelta;
            rBeam = rBeamA + w .* rBeamDelta;
            % 优化：均匀网格用等比数列相位（tau 仍传入供 stent 使用）
            if isUniform_g
                nR_loc = numel(idx);
                phaseStart = exp(-1i * omega * tau(1));
                if nR_loc > 1
                    phaseStep = exp(-1i * omega * (dr_g / drab * tauDelta));
                    phase_g = cumprod([phaseStart, repmat(phaseStep, 1, nR_loc - 1)]);
                else
                    phase_g = phaseStart;
                end
            else
                phase_g = exp(-1i * omega * tau);
            end
            lastRoot = target.lastRoots{iz};
            [linePressure, lastRoot] = local_accumulate_geometric_gaussian_line( ...
                target.pressure(iz, idx), meta, beamInfo, ...
                n, Q, tau, phase_g, c, rBeam, lastRoot, G);
            target.pressure(iz, idx) = linePressure;
            target.lastRoots{iz} = lastRoot;
        end

    otherwise
        error('Unsupported beam family code: %d (%s)', meta.familyCode, meta.family);
end
end


function [linePressure, lastRoot] = local_accumulate_geometric_hat_line( ...
    linePressure, meta, beamInfo, n, Q, phase, c, rBeam, lastRoot, G)
% phase 由调用方预计算（均匀网格用等比数列，非均匀用 exp 向量）

qReal = real(Q);
W = meta.widthScale * abs(qReal) * meta.absDeltaAlpha ./ max(abs(c), eps);
valid = isfinite(rBeam) & (rBeam > 0) & isfinite(c) & isfinite(Q) & isfinite(W) & (W > 0) ...
    & (abs(n) <= W);

if ~any(valid)
    return;
end

qEff = qReal;
qEff(qEff == 0) = eps;
shape = (W - abs(n)) ./ W;

% 优化：向量化 continuous_sqrt（实数 rootArg 特化）
% rootArg = c/(rBeam*qEff) 为实数：qEff>0 时根为实数，qEff<0 时根为纯虚数。
% signFactor (±1) 由 lastRoot 决定，同段内 qEff 符号通常不变（无焦散穿越）。
rootMag = sqrt(abs(c) ./ (rBeam .* abs(qEff)));  % 振幅模，实数，向量化
natPhase = ones(size(qEff));
natPhase(qEff < 0) = 1i;                          % qEff<0 时自然分支为纯虚数

firstIdx = find(valid, 1, 'first');
posRoot0 = rootMag(firstIdx) * natPhase(firstIdx);
if isempty(lastRoot) || abs(posRoot0 - lastRoot) <= abs(-posRoot0 - lastRoot)
    signFactor = 1;
else
    signFactor = -1;
end

rootVals = (signFactor * rootMag) .* natPhase;    % 向量化根值

linePressure(valid) = linePressure(valid) + ...
    G(valid) .* beamInfo.hatAmpScale .* rootVals(valid) .* shape(valid) .* phase(valid);

lastRoot = rootVals(find(valid, 1, 'last'));
end


function [linePressure, lastRoot] = local_accumulate_geometric_gaussian_line( ...
    linePressure, meta, beamInfo, n, Q, tau, phase, c, rBeam, lastRoot, G)
% phase 由调用方预计算；tau 仍用于 stent 宽度下限计算

qReal = real(Q);
qAbs = abs(qReal);
Wgeom = meta.widthScale * qAbs * meta.absDeltaAlpha ./ max(abs(c), eps);

if meta.useStent
    lambda = abs(c) / max(meta.freq, eps);
    Wfloor = min(0.2 * meta.freq * max(real(tau), 0.0), pi * lambda);
else
    Wfloor = zeros(size(Wgeom));
end

W = max(Wgeom, Wfloor);
valid = isfinite(rBeam) & (rBeam > 0) & isfinite(c) & isfinite(Q) & isfinite(W) & (W > 0);
if meta.autoWindow
    valid = valid & (abs(n) <= meta.gaussianWindow .* W);
end

if ~any(valid)
    return;
end

qEff = qReal;
maskFloor = W > Wgeom;
if any(maskFloor)
    qFloor = max(abs(c(maskFloor)) .* W(maskFloor) / meta.widthDenom, eps);
    qEff(maskFloor) = qFloor;
    negMask = maskFloor & (qReal < 0);
    qEff(negMask) = -max(qAbs(negMask), abs(qEff(negMask)));
    posMask = maskFloor & ~negMask;
    qEff(posMask) = max(qAbs(posMask), qEff(posMask));
end
zeroMask = (qEff == 0);
qEff(zeroMask) = eps;

A = abs(beamInfo.q0Source ./ qEff);
valid = valid & isfinite(A) & (A > 0);
if ~any(valid)
    return;
end

shape = exp(-0.5 * (n ./ W).^2) ./ (W .* A);

% 优化：向量化 continuous_sqrt（实数 rootArg 特化，同 hat 分支）
rootMag = sqrt(abs(c) ./ (rBeam .* abs(qEff)));  % 振幅模，实数，向量化
natPhase = ones(size(qEff));
natPhase(qEff < 0) = 1i;                          % qEff<0 时自然分支为纯虚数

firstIdx = find(valid, 1, 'first');
posRoot0 = rootMag(firstIdx) * natPhase(firstIdx);
if isempty(lastRoot) || abs(posRoot0 - lastRoot) <= abs(-posRoot0 - lastRoot)
    signFactor = 1;
else
    signFactor = -1;
end

rootVals = (signFactor * rootMag) .* natPhase;    % 向量化根值

linePressure(valid) = linePressure(valid) + ...
    G(valid) .* beamInfo.gaussianAmpScale .* rootVals(valid) .* shape(valid) .* phase(valid);

lastRoot = rootVals(find(valid, 1, 'last'));
end


function [linePressure, lastRoot] = local_accumulate_paraxial_line( ...
    linePressure, meta, config, beamInfo, n, P, Q, tau, c, rBeam, lastRoot, G)
% 优化：paraxial 分支向量化，替代原来对每个接收点逐一调用 local_beam_contribution_fast 的循环。
%
% rootArg = c/(rBeam*Q) 为复数。先用 MATLAB 自然分支批量计算，
% 再对 valid 子集通过向量化内积检测分支跳变并修正，无显式逐元素循环。

% 优化：缓存 omega，消除每次调用的结构体字段访问
persistent cachedOmega
if isempty(cachedOmega)
    cachedOmega = config.omega;
end
omega = cachedOmega;

pq = P ./ Q;
valid = isfinite(rBeam) & (rBeam > 0) & isfinite(c) & isfinite(pq) & (imag(pq) < 0);
if ~any(valid)
    return;
end

L = sqrt(-2.0 ./ (omega .* imag(pq)));
valid = valid & isfinite(L);
if meta.autoWindow
    valid = valid & (abs(n) <= meta.windowRadii .* L);
end
if ~any(valid)
    return;
end

% 向量化 continuous_sqrt（复数 rootArg）
% 1. 按 MATLAB 自然分支批量计算
rootVecNat = sqrt(c ./ (rBeam .* Q));

% 2. 从 lastRoot 确定初始 signFactor
validIdx = find(valid);
posRoot0  = rootVecNat(validIdx(1));
if isempty(lastRoot) || lastRoot == 0 || ...
        abs(posRoot0 - lastRoot) <= abs(-posRoot0 - lastRoot)
    signFactor = 1;
else
    signFactor = -1;
end

% 3. 向量化检测 valid 子集内的分支跳变
%    条件：Re(r_k * conj(r_{k-1})) < 0 表示需在 k 处翻转
if numel(validIdx) > 1
    vr = rootVecNat(validIdx);
    flipEvents  = real(vr(2:end) .* conj(vr(1:end-1))) < 0;
    cumFlips    = mod(cumsum([false; flipEvents(:)]), 2);
    vr(logical(cumFlips)) = -vr(logical(cumFlips));
    rootVecNat(validIdx) = vr;
end

if signFactor == -1
    rootVecNat(validIdx) = -rootVecNat(validIdx);
end

% 4. 向量化计算相位与压力贡献
phase = exp(-1i * omega .* (tau + 0.5 .* pq .* n.^2));
contribs = G .* meta.deltaAlpha .* beamInfo.Aconst .* rootVecNat .* phase;
linePressure(valid) = linePressure(valid) + contribs(valid);

lastRoot = rootVecNat(validIdx(end));
end


function [i1, i2] = local_receiver_index_span(receiver_ranges, gridInfo, rMin, rMax)
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

    i1 = max(1, ceil (x1 - 1.0e-12));
    i2 = min(nR, floor(x2 + 1.0e-12));

    if i2 < i1
        i1 = 1;
        i2 = 0;
    end
else
    i1 = find(receiver_ranges >= rMin, 1, 'first');
    i2 = find(receiver_ranges <= rMax, 1, 'last');

    if isempty(i1) || isempty(i2) || i2 < i1
        i1 = 1;
        i2 = 0;
    end
end
end


function [contrib, lastRoot, used] = local_beam_contribution_fast( ...
    meta, config, beamInfo, n, P, Q, tau, c, rBeam, lastRoot, G)

contrib = 0;
used = false;

if ~isfinite(rBeam) || rBeam <= 0 || ~isfinite(c) || ~isfinite(P) || ~isfinite(Q)
    return;
end

switch meta.familyCode
    case 1
        pq = P / Q;
        if ~(imag(pq) < 0)
            return;
        end

        L = sqrt(-2.0 / (config.omega * imag(pq)));
        if ~isfinite(L)
            return;
        end

        if meta.autoWindow && abs(n) > meta.windowRadii * L
            return;
        end

        [rootVal, lastRoot] = utils.continuous_sqrt(c / (rBeam * Q), lastRoot);
        phase   = exp(-1i * config.omega * (tau + 0.5 * pq * n^2));
        contrib = G * meta.deltaAlpha * beamInfo.Aconst * rootVal * phase;
        used    = true;

    case 2
        qReal = real(Q);
        qAbs = abs(qReal);
        W = meta.widthScale * qAbs * meta.absDeltaAlpha / max(abs(c), eps);
        if ~isfinite(W) || W <= 0 || abs(n) > W
            return;
        end

        if qReal == 0
            qEff = eps;
        else
            qEff = qReal;
        end
        shape = (W - abs(n)) / W;
        [rootVal, lastRoot] = utils.continuous_sqrt(c / (rBeam * complex(qEff, 0.0)), lastRoot);
        amp     = beamInfo.hatAmpScale * rootVal;
        phase   = exp(-1i * config.omega * tau);
        contrib = G * amp * shape * phase;
        used    = true;

    case 3
        qReal = real(Q);
        qAbs = abs(qReal);
        Wgeom = meta.widthScale * qAbs * meta.absDeltaAlpha / max(abs(c), eps);
        if meta.useStent
            lambda = abs(c) / max(meta.freq, eps);
            Wfloor = min(0.2 * meta.freq * max(real(tau), 0.0), pi * lambda);
        else
            Wfloor = 0.0;
        end
        W = max(Wgeom, Wfloor);
        if ~isfinite(W) || W <= 0
            return;
        end

        if meta.autoWindow && abs(n) > meta.gaussianWindow * W
            return;
        end

        if W > Wgeom
            qFloor = max(abs(c) * W / meta.widthDenom, eps);
            if qReal < 0
                qEff = -max(qAbs, qFloor);
            else
                qEff = max(qAbs, qFloor);
            end
        else
            if qReal == 0
                qEff = eps;
            else
                qEff = qReal;
            end
        end

        A = abs(beamInfo.q0Source / qEff);
        if ~(isfinite(A) && A > 0)
            return;
        end

        shape = exp(-0.5 * (n / W)^2) / (W * A);
        [rootVal, lastRoot] = utils.continuous_sqrt(c / (rBeam * complex(qEff, 0.0)), lastRoot);
        amp     = beamInfo.gaussianAmpScale * rootVal;
        phase   = exp(-1i * config.omega * tau);
        contrib = G * amp * shape * phase;
        used    = true;

    otherwise
        error('Unsupported beam family: %s', meta.family);
end
end
