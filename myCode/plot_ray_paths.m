function fig = plot_ray_paths(config, beams)
%PLOT_RAY_PATHS Plot 2D ray paths.

fig = figure('Color', 'w');
hold on;
for ib = 1:numel(beams)
    plot(beams{ib}.r, beams{ib}.z, 'k-', 'LineWidth', 0.9);
end
xlabel(sprintf('Range (%s)', config.units.length_name));
ylabel('Depth');
title(sprintf('%s ray trace [%s]', config.title, config.beam.family), 'Interpreter', 'none');
set(gca, 'YDir', 'reverse', 'Box', 'on');
grid on;
xlim([0, config.ray.rmax]);
ylim([0, config.ray.zmax]);
end
