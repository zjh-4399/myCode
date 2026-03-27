function hit = local_find_bottom_hit(y, yTry, ds_left, config)
%LOCAL_FIND_BOTTOM_HIT Find the earliest intersection of one beam step
%with the piecewise-linear bottom boundary.
%
% INPUT
%   y      : current beam state at the start of the substep
%            y = [r; z; pr; pz; Re(P); Im(P); Re(Q); Im(Q); tau]
%   yTry   : trial beam state after stepping ds_left
%   ds_left: current remaining arc-length step
%   config : solver config, with config.bottom preprocessed by
%            preprocess_bottom_2d(config)
%
% OUTPUT
%   hit.exists : true if the segment y -> yTry intersects the bottom
%   hit.rHit   : hit range
%   hit.zHit   : hit depth
%   hit.segId  : bottom facet index that is hit
%   hit.dsHit  : traveled arc-length from y to hit point
%   hit.nBdry  : boundary normal at hit point (filled by attach_bottom_normal)
%   hit.tBdry  : boundary tangent at hit point (filled by attach_bottom_normal)
%   hit.kappa  : boundary curvature at hit point (filled by attach_bottom_normal)

    % 优化：用 persistent 模板替代每次现场构造命中结构体，消除重复的字段分配开销
    % 切换 config 时可调用 clear local_find_bottom_hit 重置缓存
    persistent bottomEnabled cachedBottom cachedNSeg emptyHit
    if isempty(bottomEnabled)
        bottomEnabled = isfield(config, 'boundary') && isfield(config.boundary, 'bottom') ...
            && config.boundary.bottom.enabled ...
            && isfield(config, 'bottom') && ~isempty(config.bottom) ...
            && isfield(config.bottom, 'r') && isfield(config.bottom, 'z');
        if bottomEnabled
            cachedBottom = config.bottom;
            cachedNSeg   = numel(config.bottom.r) - 1;
            if cachedNSeg < 1
                bottomEnabled = false;
                cachedBottom  = [];
                cachedNSeg    = 0;
            end
        else
            cachedBottom = [];
            cachedNSeg   = 0;
        end
        emptyHit = struct( ...
            'exists', false, ...
            'rHit',   nan, ...
            'zHit',   nan, ...
            'segId',  nan, ...
            'dsHit',  inf, ...
            'nBdry',  [], ...
            'tBdry',  [], ...
            'kappa',  [] );
    end
    hit = emptyHit;
    if ~bottomEnabled
        return;
    end

    rA = y(1);
    zA = y(2);
    rB = yTry(1);
    zB = yTry(2);

    if ~all(isfinite([rA, zA, rB, zB, ds_left]))
        return;
    end

    vr = rB - rA;
    vz = zB - zA;

    % 零长度步，直接返回
    if hypot(vr, vz) <= eps
        return;
    end

    bottom = cachedBottom;
    nSeg = cachedNSeg;

    % -------------------------------------------------------------
    % Search the earliest intersection with the bottom polyline
    % -------------------------------------------------------------
    tBest = inf;
    segBest = nan;
    rBest = nan;
    zBest = nan;

    % 数值容差
    tolParallel = 1e-14;
    tolParam    = 1e-10;
    tolBox      = 1e-12;

    % 射线当前小步的包围盒
    ray_rmin = min(rA, rB);
    ray_rmax = max(rA, rB);
    ray_zmin = min(zA, zB);
    ray_zmax = max(zA, zB);

    segStart = find(bottom.r(2:end) >= ray_rmin - tolBox, 1, 'first');
    segEnd = find(bottom.r(1:end-1) <= ray_rmax + tolBox, 1, 'last');
    if isempty(segStart) || isempty(segEnd) || segEnd < segStart
        return;
    end

    % 优化：向量化底部分段交叉检测，消除候选分段上的标量循环
    segRange = segStart:segEnd;
    r1v = bottom.r(segRange);
    z1v = bottom.z(segRange);
    r2v = bottom.r(segRange + 1);
    z2v = bottom.z(segRange + 1);

    seg_rminv = min(r1v, r2v);
    seg_rmaxv = max(r1v, r2v);
    seg_zminv = min(z1v, z2v);
    seg_zmaxv = max(z1v, z2v);

    bboxOk = ~((ray_rmax < seg_rminv - tolBox) | (ray_rmin > seg_rmaxv + tolBox) | ...
                (ray_zmax < seg_zminv - tolBox) | (ray_zmin > seg_zmaxv + tolBox));

    if any(bboxOk)
        r1c = r1v(bboxOk);  z1c = z1v(bboxOk);
        wrv = r2v(bboxOk) - r1c;  wzv = z2v(bboxOk) - z1c;
        rhs_rv = r1c - rA;  rhs_zv = z1c - zA;
        denomv = vr * wzv - vz * wrv;
        nonPar = abs(denomv) >= tolParallel;
        tv = (rhs_rv .* wzv - rhs_zv .* wrv) ./ denomv;
        uv = (rhs_rv .* vz  - rhs_zv .* vr)  ./ denomv;
        hitMask = nonPar & (tv > tolParam) & (tv <= 1.0 + tolParam) & ...
                  (uv >= -tolParam) & (uv <= 1.0 + tolParam);
        if any(hitMask)
            tv(~hitMask) = Inf;
            [tRaw, bestLocal] = min(tv);
            tBest = min(max(tRaw, 0.0), 1.0);
            segCands = segRange(bboxOk);
            segBest = segCands(bestLocal);
            rBest = rA + tBest * vr;
            zBest = zA + tBest * vz;
        end
    end

    if ~isfinite(tBest)
        return;
    end

    % -------------------------------------------------------------
    % Assemble hit result
    % -------------------------------------------------------------
    hit.exists = true;
    hit.rHit   = rBest;
    hit.zHit   = zBest;
    hit.segId  = segBest;
    hit.dsHit  = ds_left * tBest;

    % 在这里补边界切向/法向/曲率信息
    hit = attach_bottom_normal(hit, config);
end
