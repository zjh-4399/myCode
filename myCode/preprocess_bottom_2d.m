function bottom = preprocess_bottom_2d(config)
%PREPROCESS_BOTTOM_2D Precompute 2D bottom geometry for reflection handling.
%
% Bellhop-like conventions used here:
%   - boundary polyline nodes are ordered by increasing range r
%   - segment tangent tSeg points in increasing-r direction
%   - bottom normal nSeg is the OUTWARD normal of the water column
%     (for a flat bottom, nSeg = [0; +1] in [r,z] coordinates)
%
% We also precompute:
%   - node tangents / normals for curvilinear-like interpolation
%   - segment curvature estimates
%
% Supported curvature estimators:
%   config.boundary.bottom.curvature_estimator =
%       'bellhop_dss'   % recommended, closest to current Bellhop code path
%       'angle'         % dphi/ds from node tangent angle
%       'zero'          % disable curvature

r = config.boundary.bottom.r(:);
z = config.boundary.bottom.z(:);

if numel(r) ~= numel(z) || numel(r) < 2
    error('boundary.bottom.r / z must have the same length >= 2.');
end
if any(diff(r) <= 0)
    error('boundary.bottom.r must be strictly increasing.');
end

nPts = numel(r);
nSegCount = nPts - 1;

dr = diff(r);
dz = diff(z);
L  = hypot(dr, dz);

if any(L <= 0)
    error('Bottom polyline contains a zero-length segment.');
end

% -------------------------------------------------------------------------
% 1) Segment tangent / normal
% -------------------------------------------------------------------------
% Tangent follows increasing range, same spirit as Bellhop Bdry(ii)%t
tSeg = [dr ./ L, dz ./ L];

% Bottom outward normal, Bellhop-style:
% for t = [tx, tz], n = [-tz, +tx]
nSeg = [-tSeg(:,2), tSeg(:,1)];

% First derivative dz/dr on each segment
DxSeg = dz ./ dr;

% -------------------------------------------------------------------------
% 2) Node tangent / normal (curvilinear-like)
% -------------------------------------------------------------------------
% Bellhop sets end-node tangents to horizontal because the bathymetry is
% effectively extended flat at the two ends.
tNode = zeros(nPts, 2);
tNode(1,:)   = [1.0, 0.0];
tNode(end,:) = [1.0, 0.0];

for i = 2:nPts-1
    v = 0.5 * tSeg(i-1,:) + 0.5 * tSeg(i,:);
    if norm(v) < 1.0e-12
        v = tSeg(i,:);
    end
    tNode(i,:) = v / norm(v);
end

nNode = [-tNode(:,2), tNode(:,1)];

% Tangent angle at nodes
phiNode = unwrap(atan2(tNode(:,2), tNode(:,1)));

% -------------------------------------------------------------------------
% 3) Curvature estimates
% -------------------------------------------------------------------------
% (a) angle-based: kappa = dphi/ds
kappaAngle = diff(phiNode) ./ max(L, eps);

% (b) Bellhop-like Dss estimator
% Bellhop currently overrides kappa by Dss = Dxx * t_x^3.
DxxSeg = zeros(nSegCount, 1);
DssSeg = zeros(nSegCount, 1);

if nSegCount >= 2
    for i = 1:nSegCount-1
        DxxSeg(i) = (DxSeg(i+1) - DxSeg(i)) / max(r(i+1) - r(i), eps);
    end
    DxxSeg(end) = DxxSeg(end-1);
    DssSeg = DxxSeg .* (tSeg(:,1) .^ 3);
end

curvEstimator = 'bellhop_dss';
if isfield(config.boundary.bottom, 'curvature_estimator') && ...
        ~isempty(config.boundary.bottom.curvature_estimator)
    curvEstimator = lower(strtrim(char(string(config.boundary.bottom.curvature_estimator))));
end

switch curvEstimator
    case 'bellhop_dss'
        kappaSeg = DssSeg;
    case 'angle'
        kappaSeg = kappaAngle;
    case 'zero'
        kappaSeg = zeros(nSegCount, 1);
    otherwise
        error('Unsupported boundary.bottom.curvature_estimator: %s', curvEstimator);
end

% -------------------------------------------------------------------------
% 4) Assemble output
% -------------------------------------------------------------------------
bottom = struct();
bottom.r = r;
bottom.z = z;

bottom.L = L;
bottom.tSeg = tSeg;
bottom.nSeg = nSeg;
bottom.DxSeg = DxSeg;

bottom.tNode = tNode;
bottom.nNode = nNode;
bottom.phiNode = phiNode;

bottom.kappaAngle = kappaAngle;
bottom.DxxSeg = DxxSeg;
bottom.DssSeg = DssSeg;
bottom.kappaSeg = kappaSeg;

bottom.interp = lower(strtrim(char(string(config.boundary.bottom.interp))));
bottom.curvature_estimator = curvEstimator;
end