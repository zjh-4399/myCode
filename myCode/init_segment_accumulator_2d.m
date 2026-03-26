function target = init_segment_accumulator_2d(config, depths, ranges, pressure, name)
%INIT_SEGMENT_ACCUMULATOR_2D Initialize one per-beam online segment accumulator.
%
% The pressure field is carried in from the global accumulator, while
% branch-tracking state (lastRoots) is reset for each newly traced beam.

if nargin < 4 || isempty(pressure)
    pressure = complex(zeros(numel(depths), numel(ranges)));
end
if nargin < 5 || isempty(name)
    name = 'target';
end

target = struct();
target.name = char(string(name));
target.depths = depths(:).';
target.ranges = ranges(:).';
target.pressure = complex(pressure);
target.lastRoots = cell(numel(target.depths), 1);
target.gridInfo = local_receiver_grid_info(target.ranges);
target.beamMeta = local_unpack_beam_meta(config);
end


function meta = local_unpack_beam_meta(config)
meta = struct();

if isfield(config, 'beam') && isfield(config.beam, 'model') && isfield(config.beam.model, 'family')
    meta.family = lower(char(string(config.beam.model.family)));
elseif isfield(config, 'beam') && isfield(config.beam, 'family')
    meta.family = lower(char(string(config.beam.family)));
else
    error('Cannot find beam family in config.beam.family or config.beam.model.family.');
end

switch meta.family
    case 'paraxial'
        meta.familyCode = 1;
    case 'geometric_hat'
        meta.familyCode = 2;
    case 'geometric_gaussian'
        meta.familyCode = 3;
    otherwise
        error('Unsupported beam family: %s', meta.family);
end

meta.autoWindow = local_read_beam_field(config, 'auto_window', false);
meta.windowRadii = local_read_beam_field(config, 'window_radii', 4.0);
meta.widthScale  = local_read_beam_field(config, 'width_scale', 1.0);
meta.widthFloor  = local_read_beam_field(config, 'width_floor', 0.0);
meta.useStent    = local_read_beam_field(config, 'use_stent', false);
meta.ibwin       = local_read_beam_field(config, 'ibwin', 4);

if isfield(config.beam, 'delta_alpha_rad') && ~isempty(config.beam.delta_alpha_rad)
    meta.deltaAlpha = config.beam.delta_alpha_rad;
else
    error('config.beam.delta_alpha_rad is missing. Make sure validate_config_2d has run.');
end
meta.absDeltaAlpha = abs(meta.deltaAlpha);
if meta.absDeltaAlpha <= 0
    meta.absDeltaAlpha = deg2rad(1.0);
end
meta.widthDenom = max(meta.widthScale * meta.absDeltaAlpha, eps);
meta.freq = config.omega / (2.0 * pi);
meta.gaussianWindow = meta.ibwin;
if isempty(meta.gaussianWindow) || ~isfinite(meta.gaussianWindow) || meta.gaussianWindow <= 0
    meta.gaussianWindow = meta.windowRadii;
end
if isempty(meta.gaussianWindow) || ~isfinite(meta.gaussianWindow) || meta.gaussianWindow <= 0
    meta.gaussianWindow = 4.0;
end
end


function value = local_read_beam_field(config, fieldName, defaultValue)
if isfield(config, 'beam') && isfield(config.beam, 'model') && isfield(config.beam.model, fieldName) ...
        && ~isempty(config.beam.model.(fieldName))
    value = config.beam.model.(fieldName);
elseif isfield(config, 'beam') && isfield(config.beam, fieldName) ...
        && ~isempty(config.beam.(fieldName))
    value = config.beam.(fieldName);
else
    value = defaultValue;
end
end


function gridInfo = local_receiver_grid_info(receiver_ranges)
gridInfo = struct();
gridInfo.isUniform = false;
gridInfo.isIncreasing = false;
gridInfo.r0 = [];
gridInfo.dr = [];
gridInfo.rMin = [];
gridInfo.rMax = [];

rr = receiver_ranges(:);
if numel(rr) <= 1
    gridInfo.isUniform = true;
    gridInfo.isIncreasing = true;
    if ~isempty(rr)
        gridInfo.r0 = rr(1);
    else
        gridInfo.r0 = 0.0;
    end
    gridInfo.dr = 1.0;
    if ~isempty(rr)
        gridInfo.rMin = rr(1);
        gridInfo.rMax = rr(end);
    else
        gridInfo.rMin = 0.0;
        gridInfo.rMax = 0.0;
    end
    return;
end

d = diff(rr);
gridInfo.isIncreasing = all(d > 0);
if ~gridInfo.isIncreasing
    return;
end

d0 = d(1);
tol = 1.0e-10 * max(1.0, abs(d0));
gridInfo.isUniform = all(abs(d - d0) <= tol);
gridInfo.r0 = rr(1);
gridInfo.dr = d0;
gridInfo.rMin = rr(1);
gridInfo.rMax = rr(end);
end
