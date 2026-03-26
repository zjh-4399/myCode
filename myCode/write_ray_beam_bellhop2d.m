function write_ray_beam_bellhop2d(fid, beam, config)
%WRITE_RAY_BEAM_BELLHOP2D Append one 2D beam to a Bellhop-style ASCII .ray file.
%
% Per-beam block:
%   alpha0
%   N2 NumTopBnc NumBotBnc
%   r z   (N2 rows)

if fid < 0
    error('Invalid file identifier for .ray writing.');
end

maxPoints = local_get_numeric(config, {'output','ray_max_points'}, 500000);
[rc, zc] = local_compress_ray_points(beam.r(:), beam.z(:), beam, config, maxPoints);

fprintf(fid, '%.16g\n', beam.alpha_deg);
fprintf(fid, '%d %d %d\n', numel(rc), local_get_bounces(beam, 'surface_bounces'), local_get_bounces(beam, 'bottom_bounces'));
for is = 1:numel(rc)
    fprintf(fid, '%.16g %.16g\n', rc(is), zc(is));
end
end


function value = local_get_numeric(config, pathParts, defaultValue)
value = defaultValue;
S = config;
for k = 1:numel(pathParts)
    fieldName = pathParts{k};
    if isstruct(S) && isfield(S, fieldName) && ~isempty(S.(fieldName))
        S = S.(fieldName);
    else
        return;
    end
end
if isnumeric(S) && isscalar(S) && isfinite(S)
    value = S;
end
end


function n = local_get_bounces(beam, fieldName)
if isfield(beam, fieldName) && ~isempty(beam.(fieldName))
    n = beam.(fieldName);
else
    n = 0;
end
end


function [rKeep, zKeep] = local_compress_ray_points(r, z, beam, config, maxPoints)
% Mirror Bellhop's WriteRay2D logic:
% keep every iSkip point, points near top/bottom, and the last point.

n = numel(r);
if n == 0
    rKeep = zeros(0,1);
    zKeep = zeros(0,1);
    return;
end
if n == 1
    rKeep = r;
    zKeep = z;
    return;
end

iSkip = max(floor(n / max(maxPoints,1)), 1);
zTop = local_top_depth(config);
zBot = local_bottom_depth(config, beam);

keep = false(n,1);
keep(1) = true;
for is = 2:n
    nearBdry = min(zBot - z(is), z(is) - zTop) < 0.2;
    if nearBdry || mod(is, iSkip) == 0 || is == n
        keep(is) = true;
    end
end

rKeep = r(keep);
zKeep = z(keep);
end


function zTop = local_top_depth(config)
zTop = 0.0;
if isfield(config, 'boundary') && isfield(config.boundary, 'surface') ...
        && isfield(config.boundary.surface, 'depth') && ~isempty(config.boundary.surface.depth)
    zTop = config.boundary.surface.depth;
end
end


function zBot = local_bottom_depth(config, beam)
if nargin >= 2 && isfield(beam, 'z') && ~isempty(beam.z)
    zBot = max(beam.z(:));
else
    zBot = 0.0;
end

if isfield(config, 'bottom') && isstruct(config.bottom)
    if isfield(config.bottom, 'z') && ~isempty(config.bottom.z)
        zBot = max(zBot, max(config.bottom.z(:)));
        return;
    end
    if isfield(config.bottom, 'zb') && ~isempty(config.bottom.zb)
        zBot = max(zBot, max(config.bottom.zb(:)));
        return;
    end
end

if isfield(config, 'boundary') && isfield(config.boundary, 'bottom')
    if isfield(config.boundary.bottom, 'z') && ~isempty(config.boundary.bottom.z)
        zBot = max(zBot, max(config.boundary.bottom.z(:)));
        return;
    end
    if isfield(config.boundary.bottom, 'depth') && ~isempty(config.boundary.bottom.depth)
        zBot = max(zBot, max(config.boundary.bottom.depth(:)));
        return;
    end
end

if zBot <= 0 && isfield(config, 'ray') && isfield(config.ray, 'zmax')
    zBot = config.ray.zmax;
end
end
