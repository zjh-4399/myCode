function results = run_solver_2d(config)
%RUN_SOLVER_2D Main 2D solver entry point with Bellhop-style config.
%   results = RUN_SOLVER_2D(config)
%
% Key behaviors in this version:
%   1) coherent_tl path uses segment-level online accumulation;
%   2) rays can be streamed to a Bellhop-style .ray file during tracing;
%   3) beam cell storage is controlled explicitly by config.output.store_beams.

config = validate_config_2d(config);

results = struct();
results.config = config;
results.output_type = config.output.model.type;

storeBeams = local_should_store_beams(config);
rayWriter = local_open_ray_writer(config);
rayCleanup = onCleanup(@() local_close_ray_writer(rayWriter)); %#ok<NASGU>

switch config.output.model.type
    case 'ray'
        traceOpts = struct('storeTrace', true);
        if storeBeams
            beams = cell(config.ray.nbeams, 1);
        else
            beams = {};
        end

        for ib = 1:config.ray.nbeams
            beam = trace_beam_auxiliary(config, config.ray.theta_list_deg(ib), [], traceOpts);
            if rayWriter.enabled
                write_ray_beam_bellhop2d(rayWriter.fid, beam, config);
            end
            if storeBeams
                beams{ib} = beam;
            end
        end

        if storeBeams
            results.beams = beams;
            results.rays = beams;
        else
            results.beams = {};
            results.rays = {};
        end

        if rayWriter.enabled
            results.ray_file = rayWriter.filepath;
            local_close_ray_writer(rayWriter);
            rayWriter.enabled = false;
        end

        if config.output.plot_ray
            results.fig_ray = local_plot_rays(config, results);
        end

    case 'coherent_tl'
        fieldEnabled = isfield(config, 'field') && isfield(config.field, 'enabled') && config.field.enabled;
        onlineTargets = local_init_online_targets(config, fieldEnabled);
        needTraceStorage = storeBeams || rayWriter.enabled;
        traceOpts = struct('storeTrace', needTraceStorage);

        if storeBeams
            beams = cell(config.ray.nbeams, 1);
        else
            beams = {};
        end

        for ib = 1:config.ray.nbeams
            onlineTargets = local_reset_online_targets(onlineTargets);

            [beam, onlineTargets] = trace_beam_auxiliary( ...
                config, config.ray.theta_list_deg(ib), onlineTargets, traceOpts);

            if rayWriter.enabled
                write_ray_beam_bellhop2d(rayWriter.fid, beam, config);
            end
            if storeBeams
                beams{ib} = beam;
            end
        end

        receiverPressure = onlineTargets(1).pressure;
        results.pressure = receiverPressure;
        results.tl = pressure_to_tl(receiverPressure);
        results.receiver_depths = config.receiver.depths;
        results.receiver_ranges = config.receiver.ranges;

        if storeBeams
            results.beams = beams;
        else
            results.beams = {};
        end

        if rayWriter.enabled
            results.ray_file = rayWriter.filepath;
            local_close_ray_writer(rayWriter);
            rayWriter.enabled = false;
        end

        if fieldEnabled
            fieldPressure = onlineTargets(2).pressure;
            results.field_pressure = fieldPressure;
            results.field_tl = pressure_to_tl(fieldPressure);
            results.field_depths = config.field.depths;
            results.field_ranges = config.field.ranges;
        end

        if config.output.plot_ray
            results.fig_ray = local_plot_rays(config, results);
        end
        if config.output.plot_tl
            results.fig_tl = plot_tl_curve(config, config.receiver.ranges, results.tl(1,:));
        end
        if isfield(results, 'field_tl') && config.output.plot_tl
            results.fig_field = plot_tl_field(config, results.field_ranges, ...
                results.field_depths, results.field_tl);
        end

    otherwise
        error('Unsupported output model type: %s', config.output.model.type);
end
end


function targets = local_init_online_targets(config, fieldEnabled)
targets = init_segment_accumulator_2d( ...
    config, config.receiver.depths, config.receiver.ranges, [], 'receiver');

if fieldEnabled
    targets(2) = init_segment_accumulator_2d( ...
        config, config.field.depths, config.field.ranges, [], 'field');
end
end


function targets = local_reset_online_targets(targets)
for it = 1:numel(targets)
    targets(it).lastRoots(:) = {[]};
end
end


function tf = local_should_store_beams(config)
% Explicit switch. When false, do not retain the beam cell array in results.
% Exception: if the user wants to plot rays immediately but does not save a
% ray file, we still need in-memory beams for plotting.

tf = local_get_output_flag(config, 'store_beams', false);
if tf
    return;
end

plotRay = local_get_output_flag(config, 'plot_ray', false);
saveRay = local_get_output_flag(config, 'save_ray', false);
if plotRay && ~saveRay
    tf = true;
end
end


function fig = local_plot_rays(config, results)
if isfield(results, 'beams') && ~isempty(results.beams)
    fig = plot_ray_paths(config, results.beams);
    return;
end

if isfield(results, 'ray_file') && ~isempty(results.ray_file)
    rayData = read_ray_file_bellhop2d(results.ray_file);
    fig = plot_ray_file_bellhop2d(rayData, config);
    return;
end

warning('plot_ray requested but neither beams nor ray_file are available; skipping ray plot.');
fig = [];
end


function writer = local_open_ray_writer(config)
writer = struct('enabled', false, 'fid', -1, 'filepath', '');
if ~local_get_output_flag(config, 'save_ray', false)
    return;
end

filepath = config.output.ray_file_path;
[fid, msg] = fopen(filepath, 'wt');
if fid < 0
    error('Failed to open ray file for writing: %s (%s)', filepath, msg);
end

write_ray_header_bellhop2d(fid, config);
writer.enabled = true;
writer.fid = fid;
writer.filepath = filepath;
end


function local_close_ray_writer(writer)
if isstruct(writer) && isfield(writer, 'enabled') && writer.enabled ...
        && isfield(writer, 'fid') && writer.fid > 0
    try
        fclose(writer.fid);
    catch
        % ignore repeated close attempts
    end
end
end


function tf = local_get_output_flag(config, fieldName, defaultValue)
if isfield(config, 'output') && isfield(config.output, fieldName) ...
        && ~isempty(config.output.(fieldName))
    tf = logical(config.output.(fieldName));
else
    tf = logical(defaultValue);
end
end
