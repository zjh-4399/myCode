function write_ray_header_bellhop2d(fid, config)
%WRITE_RAY_HEADER_BELLHOP2D Write a Bellhop-style ASCII .ray header.
%
% Format aligned to Bellhop's formatted .ray output:
%   'title'
%   freq
%   NSx NSy NSz
%   Nalpha Nbeta
%   top_depth
%   bottom_depth
%   'rz'

titleStr = local_trim_title(config);
fprintf(fid, '''%s''\n', titleStr);
fprintf(fid, '%.16g\n', config.source.freq_hz);
fprintf(fid, '%d %d %d\n', 1, 1, 1);
fprintf(fid, '%d %d\n', config.ray.nbeams, 1);
fprintf(fid, '%.16g\n', local_top_depth(config));
fprintf(fid, '%.16g\n', local_bottom_depth(config));
fprintf(fid, '''rz''\n');
end


function titleStr = local_trim_title(config)
if isfield(config, 'title') && ~isempty(config.title)
    titleStr = char(config.title);
else
    titleStr = '2D beam tracing case';
end
if numel(titleStr) > 50
    titleStr = titleStr(1:50);
end
end


function zTop = local_top_depth(config)
zTop = 0.0;
if isfield(config, 'boundary') && isfield(config.boundary, 'surface') ...
        && isfield(config.boundary.surface, 'depth') && ~isempty(config.boundary.surface.depth)
    zTop = config.boundary.surface.depth;
end
end


function zBot = local_bottom_depth(config)
if isfield(config, 'bottom') && isstruct(config.bottom)
    if isfield(config.bottom, 'z') && ~isempty(config.bottom.z)
        zBot = max(config.bottom.z(:));
        return;
    end
    if isfield(config.bottom, 'zb') && ~isempty(config.bottom.zb)
        zBot = max(config.bottom.zb(:));
        return;
    end
end

if isfield(config, 'boundary') && isfield(config.boundary, 'bottom')
    if isfield(config.boundary.bottom, 'z') && ~isempty(config.boundary.bottom.z)
        zBot = max(config.boundary.bottom.z(:));
        return;
    end
    if isfield(config.boundary.bottom, 'depth') && ~isempty(config.boundary.bottom.depth)
        zBot = max(config.boundary.bottom.depth(:));
        return;
    end
end

if isfield(config, 'ray') && isfield(config.ray, 'zmax') && ~isempty(config.ray.zmax)
    zBot = config.ray.zmax;
else
    zBot = 0.0;
end
end
