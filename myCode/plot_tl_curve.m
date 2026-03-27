function fig = plot_tl_curve(config, ranges, tlCurve)
%PLOT_TL_CURVE Plot a single TL curve.

fig = figure('Color', 'w');
plot(ranges, tlCurve, 'k-', 'LineWidth', 1.2);
grid on;
box on;
set(gca, 'YDir', 'reverse');
xlabel(sprintf('Range (%s)', config.units.length_name));
ylabel('Transmission loss (dB)');
title(sprintf('%s TL [%s]', config.title, config.beam.model.family), 'Interpreter', 'none');
xlim([min(ranges), max(ranges)]);
end
