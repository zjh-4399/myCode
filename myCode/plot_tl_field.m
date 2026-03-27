function fig = plot_tl_field(config, ranges, depths, tlField)
%PLOT_TL_FIELD Plot a 2D TL field.

fig = figure('Color', 'w');
imagesc(ranges, depths, tlField);
axis xy;
set(gca, 'YDir', 'reverse');
xlabel(sprintf('Range (%s)', config.units.length_name));
ylabel('Depth');
title(sprintf('%s TL field [%s]', config.title, config.beam.model.family), 'Interpreter', 'none');
colorbar;
colormap(flipud(jet))
caxis([60 140])
box on;
end
