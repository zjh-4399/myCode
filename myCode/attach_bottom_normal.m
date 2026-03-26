function hit = attach_bottom_normal(hit, config)
%ATTACH_BOTTOM_NORMAL Attach boundary tangent / normal / curvature at hit point.
%
% Conventions:
%   tBdry : boundary tangent
%   nBdry : OUTWARD boundary normal for the bottom
%   kappa : boundary curvature used in the reflection jump formula
%
% Supported interpolation modes:
%   'facet'         : use piecewise-linear facet tangent/normal
%   'normal_linear' : interpolate node tangents, then reconstruct normal
%
% Notes:
%   - For 'facet', default kappa = 0 (piecewise-linear facet has zero
%     curvature except at corners).
%   - For 'normal_linear', use precomputed segment curvature kappaSeg(seg).

bottom = config.bottom;

seg = hit.segId;
r1 = bottom.r(seg);
r2 = bottom.r(seg+1);

u = (hit.rHit - r1) / max(r2 - r1, eps);
u = min(max(u, 0.0), 1.0);

switch lower(bottom.interp)
    case 'facet'
        tB = bottom.tSeg(seg,:).';
        nB = bottom.nSeg(seg,:).';

        useFacetCurv = false;
        if isfield(config.boundary.bottom, 'use_curvature_on_facet') && ...
                ~isempty(config.boundary.bottom.use_curvature_on_facet)
            useFacetCurv = logical(config.boundary.bottom.use_curvature_on_facet);
        end

        if useFacetCurv
            kappa = bottom.kappaSeg(seg);
        else
            kappa = 0.0;
        end

    case 'normal_linear'
        % Bellhop curvilinear spirit:
        % interpolate node tangents, normalize, then reconstruct outward normal
        tTmp = (1.0 - u) * bottom.tNode(seg,:) + u * bottom.tNode(seg+1,:);
        if norm(tTmp) < 1.0e-12
            tTmp = bottom.tSeg(seg,:);
        end

        tB = tTmp(:) / max(norm(tTmp), eps);
        nB = [-tB(2); tB(1)];   % outward normal for bottom
        nB = nB / max(norm(nB), eps);

        kappa = bottom.kappaSeg(seg);

    otherwise
        error('Unsupported bottom.interp: %s', bottom.interp);
end

hit.tBdry = tB;
hit.nBdry = nB;
hit.kappa = kappa;
hit.uBdry = u;
end