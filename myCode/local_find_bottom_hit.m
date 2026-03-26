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

    hit = struct( ...
        'exists', false, ...
        'rHit',   nan, ...
        'zHit',   nan, ...
        'segId',  nan, ...
        'dsHit',  inf, ...
        'nBdry',  [], ...
        'tBdry',  [], ...
        'kappa',  [] );

    % -------------------------------------------------------------
    % Basic guards
    % -------------------------------------------------------------
    if ~isfield(config, 'boundary') || ~isfield(config.boundary, 'bottom') || ...
            ~config.boundary.bottom.enabled
        return;
    end

    if ~isfield(config, 'bottom') || isempty(config.bottom) || ...
            ~isfield(config.bottom, 'r') || ~isfield(config.bottom, 'z')
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

    bottom = config.bottom;
    nSeg = numel(bottom.r) - 1;
    if nSeg < 1
        return;
    end

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

    for segId = segStart:segEnd
        r1 = bottom.r(segId);
        z1 = bottom.z(segId);
        r2 = bottom.r(segId + 1);
        z2 = bottom.z(segId + 1);

        % -------------------------
        % Fast bbox rejection
        % -------------------------
        seg_rmin = min(r1, r2);
        seg_rmax = max(r1, r2);
        seg_zmin = min(z1, z2);
        seg_zmax = max(z1, z2);

        if (ray_rmax < seg_rmin - tolBox) || (ray_rmin > seg_rmax + tolBox) || ...
           (ray_zmax < seg_zmin - tolBox) || (ray_zmin > seg_zmax + tolBox)
            continue;
        end

        % -------------------------
        % Segment-segment intersection
        %
        % Ray step    : A + t * v,  t in [0,1]
        % Bottom facet: C + u * w,  u in [0,1]
        % -------------------------
        wr = r2 - r1;
        wz = z2 - z1;

        rhs_r = r1 - rA;
        rhs_z = z1 - zA;

        denom = vr * wz - vz * wr;

        % 平行/近似平行：先跳过
        % （初版足够；后面若要处理共线特殊情形再补）
        if abs(denom) < tolParallel
            continue;
        end

        t = (rhs_r * wz - rhs_z * wr) / denom;
        u = (rhs_r * vz - rhs_z * vr) / denom;

        % 命中判据：
        % 1) t 要在当前小步内部
        % 2) u 要在该 bottom facet 内
        %
        % 这里要求 t > tolParam 而不是 t >= 0，
        % 是为了避免刚在边界反射后，又因为数值舍入立刻“再次撞到同一点”
        if (t > tolParam) && (t <= 1.0 + tolParam) && ...
           (u >= -tolParam) && (u <= 1.0 + tolParam)

            if t < tBest
                tClamped = min(max(t, 0.0), 1.0);
                tBest = tClamped;
                segBest = segId;
                rBest = rA + tClamped * vr;
                zBest = zA + tClamped * vz;
            end
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
