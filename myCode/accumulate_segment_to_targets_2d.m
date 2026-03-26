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

rMinAll = min(raAll, rbAll);
rMaxAll = max(raAll, rbAll);
maskOverlap = maskFinite & (rMaxAll >= rr(1)) & (rMinAll <= rr(end));
if ~any(maskOverlap)
    return;
end

depthIdx = find(maskOverlap);
for ii = 1:numel(depthIdx)
    iz = depthIdx(ii);
    ra = raAll(iz);
    rb = rbAll(iz);
    na = naAll(iz);
    nb = nbAll(iz);

    if ~(isfinite(ra) && isfinite(rb) && isfinite(na) && isfinite(nb))
        continue;
    end

    drab = rb - ra;
    if abs(drab) < 1.0e-12
        continue;
    end

    rMin = min(ra, rb);
    rMax = max(ra, rb);
    [i1, i2] = local_receiver_index_span(rr, target.gridInfo, rMin, rMax);
    if i2 < i1
        continue;
    end

    idx = i1:i2;
    rLine = rr(idx);
    w = (rLine - ra) ./ drab;

    G     = (1.0 - w) .* pointA.gain + w .* pointB.gain;
    n     = (1.0 - w) .* na + w .* nb;
    P     = (1.0 - w) .* pointA.P + w .* pointB.P;
    Q     = (1.0 - w) .* pointA.Q + w .* pointB.Q;
    tau   = (1.0 - w) .* pointA.tau + w .* pointB.tau;
    c     = (1.0 - w) .* pointA.c + w .* pointB.c;
    rBeam = (1.0 - w) .* pointA.r + w .* pointB.r;

    lastRoot = target.lastRoots{iz};
    switch meta.familyCode
        case 2
            [linePressure, lastRoot] = local_accumulate_geometric_hat_line( ...
                target.pressure(iz, idx), meta, config, beamInfo, ...
                n, Q, tau, c, rBeam, lastRoot, G);
            target.pressure(iz, idx) = linePressure;

        case 3
            [linePressure, lastRoot] = local_accumulate_geometric_gaussian_line( ...
                target.pressure(iz, idx), meta, config, beamInfo, ...
                n, Q, tau, c, rBeam, lastRoot, G);
            target.pressure(iz, idx) = linePressure;

        otherwise
            for k = 1:numel(idx)
                [contrib, lastRoot, used] = local_beam_contribution_fast( ...
                    meta, config, beamInfo, ...
                    n(k), P(k), Q(k), tau(k), c(k), rBeam(k), lastRoot, G(k));

                if used
                    target.pressure(iz, idx(k)) = target.pressure(iz, idx(k)) + contrib;
                end
            end
    end
    target.lastRoots{iz} = lastRoot;
end
end


function [linePressure, lastRoot] = local_accumulate_geometric_hat_line( ...
    linePressure, meta, config, beamInfo, n, Q, tau, c, rBeam, lastRoot, G)

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
phase = exp(-1i * config.omega * tau);
rootArg = c ./ (rBeam .* complex(qEff, 0.0));

for k = 1:numel(linePressure)
    if ~valid(k)
        continue;
    end

    rootVal = sqrt(rootArg(k));
    if ~isempty(lastRoot) && abs(rootVal - lastRoot) > abs(-rootVal - lastRoot)
        rootVal = -rootVal;
    end
    lastRoot = rootVal;

    linePressure(k) = linePressure(k) + ...
        G(k) * beamInfo.hatAmpScale * rootVal * shape(k) * phase(k);
end
end


function [linePressure, lastRoot] = local_accumulate_geometric_gaussian_line( ...
    linePressure, meta, config, beamInfo, n, Q, tau, c, rBeam, lastRoot, G)

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
phase = exp(-1i * config.omega * tau);
rootArg = c ./ (rBeam .* complex(qEff, 0.0));

for k = 1:numel(linePressure)
    if ~valid(k)
        continue;
    end

    rootVal = sqrt(rootArg(k));
    if ~isempty(lastRoot) && abs(rootVal - lastRoot) > abs(-rootVal - lastRoot)
        rootVal = -rootVal;
    end
    lastRoot = rootVal;

    linePressure(k) = linePressure(k) + ...
        G(k) * beamInfo.gaussianAmpScale * rootVal * shape(k) * phase(k);
end
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
