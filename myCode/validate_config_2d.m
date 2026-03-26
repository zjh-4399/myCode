function config = validate_config_2d(config)
%VALIDATE_CONFIG_2D Validate a Bellhop-style 2D config and create solver aliases.
% 支持 V4 版本的 tabular 离散剖面计算及边界自适应检查

if ~isfield(config, 'dimension') || ~strcmpi(config.dimension, '2D')
    error('This project only supports 2D configurations.');
end

if ~isfield(config, 'bellhop')
    config.bellhop = struct();
end
config = apply_bellhop_option_strings_2d(config);

if ~isfield(config.beam, 'use_stent')
    config.beam.use_stent = true;
end

config.beam.model = parse_beam_options_2d(config);
config.output.model = parse_output_options_2d(config);


if ~isfield(config.output, 'store_beams') || isempty(config.output.store_beams)
    config.output.store_beams = false;
end
if ~isfield(config.output, 'ray_max_points') || isempty(config.output.ray_max_points)
    config.output.ray_max_points = 500000;
end
if ~isfield(config.output, 'ray_file_name') || isempty(config.output.ray_file_name)
    config.output.ray_file_name = [lower(strrep(config.title, ' ', '_')) '.ray'];
end

if strcmpi(config.receiver.grid_type, 'irregular')
    if ~isvector(config.receiver.depths) || ~isvector(config.receiver.ranges)
        error('For grid_type = irregular, receiver.depths and receiver.ranges must be vectors.');
    end
elseif ~strcmpi(config.receiver.grid_type, 'rectilinear')
    error('Unsupported receiver.grid_type: %s', config.receiver.grid_type);
end

% =========================================================================
% 1. V4 核心特性：处理 Tabular 离散声速剖面的预计算 (三次样条插值)
% =========================================================================
if isfield(config.env.profile, 'type') && strcmpi(config.env.profile.type, 'tabular')
    if ~isfield(config.env.profile, 'zw') || ~isfield(config.env.profile, 'cw')
        error('For tabular profile type, config.env.profile.zw and .cw must be provided.');
    end
    zw = config.env.profile.zw(:);
    cw = config.env.profile.cw(:);
    
    % 确保深度单调递增 (插值函数的严格要求)
    if any(diff(zw) <= 0)
        error('Tabular profile depths (zw) must be strictly increasing.');
    end
    
    % 构建三次样条插值多项式 (Cubic Spline)
    ppc = spline(zw, cw);
    
    % 使用无依赖的 unmkpp/mkpp 提取系数并手动求导
    [breaks, coefs, l, k, ~] = unmkpp(ppc);
    
    % 计算一阶导数系数
    coefs_z = coefs(:, 1:k-1) .* repmat((k-1:-1:1), l, 1);
    ppcz = mkpp(breaks, coefs_z);
    
    % 计算二阶导数系数
    if k > 2
        coefs_zz = coefs_z(:, 1:k-2) .* repmat((k-2:-1:1), l, 1);
        ppczz = mkpp(breaks, coefs_zz);
    else
        ppczz = mkpp(breaks, zeros(l, 1)); % 线性情况下的二阶导数为 0
    end
    
    % 存入 params 供运行时快速调用，避免 ODE 循环内重复耗时
    config.env.profile.params.ppc = ppc;
    config.env.profile.params.ppcz = ppcz;
    config.env.profile.params.ppczz = ppczz;

    if isfield(config.env.profile, 'tabular_mode') && ...
            strcmpi(char(string(config.env.profile.tabular_mode)), 'piecewise_linear')
        config.env.profile.params.pwl_zw = zw;
        config.env.profile.params.pwl_cw = cw;
        config.env.profile.params.pwl_slope = diff(cw) ./ diff(zw);
        config.env.profile.params.pwl_nseg = numel(zw) - 1;
    end
end

if isfield(config, 'numerics') && isfield(config.numerics, 'ssp_interface_jump') ...
        && config.numerics.ssp_interface_jump ...
        && isfield(config, 'env') && isfield(config.env, 'profile') ...
        && isfield(config.env.profile, 'type') && strcmpi(config.env.profile.type, 'tabular') ...
        && isfield(config.env.profile, 'tabular_mode') ...
        && ~strcmpi(char(string(config.env.profile.tabular_mode)), 'piecewise_linear')
    warning(['ssp_interface_jump only matches Bellhop weak-interface logic for ', ...
        'tabular_mode = piecewise_linear; disabling it for the current tabular_mode.']);
    config.numerics.ssp_interface_jump = false;
end

% =========================================================================
% 2. 动态修正追踪范围界限 
% 避免 generate_config 中由于 Name-Value 覆盖导致的范围截断不匹配
% =========================================================================
max_range_needed = max(config.receiver.ranges);
max_depth_needed = max(config.receiver.depths);
if isfield(config, 'field') && config.field.enabled
    max_range_needed = max([max_range_needed, max(config.field.ranges)]);
    max_depth_needed = max([max_depth_needed, max(config.field.depths)]);
end

% 如果初始设置的 rmax 或 zmax 小于所需的网格最大范围，自动扩充以覆盖计算区域 (留 5% 余量)
config.ray.rmax = max(config.ray.rmax, max_range_needed * 1.05);
config.ray.zmax = max(config.ray.zmax, max_depth_needed * 1.05);

% 如果有海深剖面，确保射线计算最深不超过环境定义的最深点
if isfield(config.env.profile, 'zw')
    config.ray.zmax = min(config.ray.zmax, max(config.env.profile.zw));
end

% 动态修正离散步长 (如果用户没有手动输入 step，且默认生成的 ds 过大或过小，在此做保护)
if isempty(config.ray.ds) || config.ray.ds <= 0
    config.ray.ds = config.ray.zmax / 25; 
end
% =========================================================================

% 射线扇面角度生成逻辑
if isempty(config.ray.theta_list_deg)
    if config.ray.auto_nbeams
        span = abs(config.ray.theta_max_deg - config.ray.theta_min_deg);
        config.ray.nbeams = max(31, ceil(6 * span));
    end
    config.ray.theta_list_deg = linspace(config.ray.theta_min_deg, ...
        config.ray.theta_max_deg, config.ray.nbeams);
else
    config.ray.theta_list_deg = config.ray.theta_list_deg(:).';
    config.ray.nbeams = numel(config.ray.theta_list_deg);
    config.ray.theta_min_deg = min(config.ray.theta_list_deg);
    config.ray.theta_max_deg = max(config.ray.theta_list_deg);
end

if numel(config.ray.theta_list_deg) > 1
    config.beam.delta_alpha_rad = deg2rad(config.ray.theta_list_deg(2) - config.ray.theta_list_deg(1));
else
    config.beam.delta_alpha_rad = deg2rad(1.0);
end

config.omega = 2*pi*config.source.freq_hz;

% Internal aliases kept for backward compatibility with the original files.
config.case_id = lower(strrep(config.title, ' ', '_'));
config.profile.type = config.env.profile.type;
config.profile = copy_profile_params(config.profile, config.env.profile.params);
config.frequency_Hz = config.source.freq_hz;

config.surface = config.boundary.surface;
config.surface.curvature_factor = config.beam.model.curvature_factor;
config.receiver.nRanges = numel(config.receiver.ranges);
config.receiver.nDepths = numel(config.receiver.depths);

config.stop.max_range = config.ray.rmax;
config.stop.max_depth = config.ray.zmax;
if ~isfield(config, 'stop') || ~isfield(config.stop, 'min_gain') || isempty(config.stop.min_gain)
    config.stop.min_gain = 5.0e-3;   % Bellhop-style amplitude cutoff
end

config.integrator.method = config.numerics.integrator;
config.integrator.ds = config.ray.ds;

% 预测积分步数
if isempty(config.numerics.max_steps)
    config.integrator.nSteps = ceil(2.0 * config.ray.rmax / max(config.ray.ds, eps)) + ...
        5 * max(1, config.boundary.surface.max_reflections);
else
    config.integrator.nSteps = config.numerics.max_steps;
end
config.integrator.total_s = config.integrator.ds * config.integrator.nSteps;

ensure_dir(config.output.results_dir);
ensure_dir(config.output.fig_dir);
config.output.ray_file_path = fullfile(config.output.results_dir, config.output.ray_file_name);

if config.boundary.bottom.enabled
    config.bottom = preprocess_bottom_2d(config);
else
    config.bottom = struct();
end
end

function profile = copy_profile_params(profile, params)
fields = fieldnames(params);
for k = 1:numel(fields)
    profile.(fields{k}) = params.(fields{k});
end
end

function ensure_dir(folderPath)
if ~exist(folderPath, 'dir')
    mkdir(folderPath);
end
end
