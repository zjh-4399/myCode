function fig = plot_ray_file_bellhop2d(rayData, config)
%PLOT_RAY_FILE_BELLHOP2D Plot rays from a Bellhop-style 2D .ray file.

if ischar(rayData) || isstring(rayData)
    rayData = read_ray_file_bellhop2d(char(rayData));
end

fig = figure('Color', 'w');
ax = axes('Parent', fig);
hold(ax, 'on');

for ib = 1:numel(rayData.beams)
    beam = rayData.beams{ib};
    plot(ax, beam.r, beam.z, 'b-');
end

set(ax, 'YDir', 'reverse');
box(ax, 'on');
grid(ax, 'on');
xlabel(ax, 'Range (m)');
ylabel(ax, 'Depth (m)');

if nargin >= 2 && isstruct(config) && isfield(config, 'title') && ~isempty(config.title)
    title(ax, sprintf('%s | rays from file', config.title), 'Interpreter', 'none');
else
    title(ax, sprintf('%s | rays from file', rayData.title), 'Interpreter', 'none');
end

if isfield(rayData, 'bottom_depth') && ~isempty(rayData.bottom_depth)
    xlim(ax, [min(local_collect_ranges(rayData.beams)) max(local_collect_ranges(rayData.beams))]);
    ylim(ax, [rayData.top_depth rayData.bottom_depth]);
end
end


function rr = local_collect_ranges(beams)
rr = [];
for k = 1:numel(beams)
    if isfield(beams{k}, 'r') && ~isempty(beams{k}.r)
        rr = [rr; beams{k}.r(:)]; %#ok<AGROW>
    end
end
if isempty(rr)
    rr = [0; 1];
end
end
