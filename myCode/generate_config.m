function config = generate_config(varargin)
%GENERATE_CONFIG针对 V4 版本的全局配置生成函数
%   自动加载所有默认参数。如果传入 Name-Value 对，则覆盖对应的默认值。
%   示例：
%   config = generate_config_v4(); % 全部使用默认值
%   config = generate_config_v4('cw', [1500, 1520], 'zw', [0, 1000]); % 覆盖 cw 和 zw
%   config = generate_config_v4('ray.nbeams', 100, 'source.depth', 500); % 覆盖嵌套属性

%% 1. 初始化所有的默认配置 (Default Settings)
% 基础信息
config.title = 'gbt_ocean_matlab_v4 default';
config.dimension = '2D';
config.reference = 'V4 Config Generator';
config.notes = '';

% 环境与声速剖面 (默认设为 tabular 离散模式)
config.env.profile.type = 'tabular';
config.env.profile.zw = [0; 5000];
config.env.profile.cw = [1500; 1500];
config.env.profile.tabular_mode = 'piecewise_linear';   % 'piecewise_linear' / 'pp_spline'

% 声源与接收器
config.source.range = 0.0;
config.source.depth = 1000.0;
config.source.freq_hz = 50.0;
config.source.beam_pattern = 'omni';

config.receiver.depths = 1000.0;
config.receiver.ranges = linspace(100, 100000.0, 10000);
config.receiver.grid_type = 'rectilinear';

% 输出模型参数
config.output.type = 'coherent_tl';
config.output.component = 'pressure';
config.output.save_ray = true;
config.output.save_field = true;
config.output.plot_ray = true;
config.output.plot_tl = true;
config.output.store_beams = false;
config.output.ray_max_points = 500000;
config.output.ray_file_name = 'beam_trace.ray';
config.output.results_dir = fullfile(pwd, 'results');
config.output.fig_dir = fullfile(pwd, 'figures');

% 声场绘制网格参数
config.field.enabled = true;
config.field.depths = 5:5:5000;
config.field.ranges = 100:100:1e5;

% 射线追踪参数
config.ray.nbeams = 88;
config.ray.auto_nbeams = false;
config.ray.theta_min_deg = -90.0;
config.ray.theta_max_deg = 90.0;
config.ray.theta_list_deg = [];
config.ray.ds = max(config.env.profile.zw)/10;
config.ray.zmax = max(config.field.depths) * 1.1;
config.ray.rmax = max(config.field.ranges) * 1.1;


% 波束参数 (Bellhop风格)
config.beam.family = 'paraxial';
config.beam.approx_code = 'C';
config.beam.init_code = 'C';
config.beam.curvature_code = 'S';
config.beam.shift_enabled = false;
config.beam.width_floor = 5.0;
config.beam.width_scale = 1.0;
config.beam.fixed_width = [];
config.beam.manual_epsilon = [];
config.beam.epmult = 1.0;
config.beam.window_radii = 5.0;
config.beam.auto_window = true;
config.beam.use_stent = true;
config.beam.rloop = 1;
config.beam.isingl = 0;
config.beam.nimage = 0;
config.beam.ibwin = 5;
config.stop.min_gain = 5.0e-3;

% 边界条件参数
config.boundary.surface.type = 'pressure_release';
config.boundary.surface.model = 'grazing_empirical';
config.boundary.surface.enabled = true;
config.boundary.surface.max_reflections = 50;

config.boundary.bottom.type = 'none';
config.boundary.bottom.model = 'none';
config.boundary.bottom.enabled = true;
config.boundary.bottom.type = 'bathymetry';

config.boundary.bottom.interp = 'normal_linear';   % facet / normal_linear
config.boundary.bottom.curvature_estimator = 'bellhop_dss';  % bellhop_dss / angle / zero
config.boundary.bottom.use_curvature_on_facet = false;

config.boundary.bottom.r = [0 config.ray.rmax];   % m
config.boundary.bottom.z = [max(config.env.profile.zw) max(config.env.profile.zw)];   % m
config.boundary.bottom.max_reflections = 50;

config.boundary.bottom.material.model = 'fluid_halfspace';
config.boundary.bottom.material.rho = 1.8;        % g/cm^3
config.boundary.bottom.material.c = 1600;         % m/s
config.boundary.bottom.material.attn = 0.8;       % dB/wavelength
config.boundary.bottom.reflection = 'plane_wave';

% 数值计算配置
config.numerics.integrator = 'rk4_fixed';   % 'rk4_fixed' / 'bellhop2_fixed'
config.numerics.branch_tracking = 'continuous_sqrt';
config.numerics.interp_scheme = 'linear';
config.numerics.max_steps = [];
config.numerics.store_all_steps = true;

config.numerics.ssp_interface_jump = true;  % Bellhop-like weak-interface jump
config.numerics.ssp_event_tol = 1e-10;

% Bellhop 选项字符串
config.bellhop.options3 = 'CCRRR';
config.bellhop.options4 = 'CS';
config.bellhop.component = '';

% 物理单位
config.units.length = 'm';
config.units.length_name = 'meter';
config.units.speed = 'm/s';


%% 2. 解析用户输入并覆盖默认值 (Override Default Settings)
if nargin > 0
    % 确保输入是成对出现的 (Name-Value)
    if mod(nargin, 2) ~= 0
        error('Input arguments must be Name-Value pairs.');
    end

    for i = 1:2:nargin
        paramName = varargin{i};
        paramValue = varargin{i+1};

        % 判断映射逻辑
        switch lower(paramName)
            % ----------------- 常用快捷映射 -----------------
            case 'cw'
                config.env.profile.cw = paramValue(:);
            case 'run_type'
                config.output.type = paramValue;
            case 'zw'
                config.env.profile.zw = paramValue(:);
            case 'freq'
                config.source.freq_hz = paramValue;
            case 'sd'        % Source Depth
                config.source.depth = paramValue;
            case 'rd'        % Receiver Depths
                config.receiver.depths = paramValue;
            case 'rr'        % Receiver Ranges
                config.receiver.ranges = paramValue(:);
            case 'fd'        % Field Depths
                config.field.depths = paramValue(:);
            case 'fr'        % Field Ranges
                config.field.ranges = paramValue(:);
            case 'nbeams'    % Number of beams
                config.ray.nbeams = paramValue(:);
            case 'theta'     % 出射角度
                config.ray.theta_min_deg = paramValue(1);
                config.ray.theta_max_deg = paramValue(2);
            case 'beam'      % 波束类型
                beamType = lower(string(paramValue));
                config.beam.family = char(beamType);

                switch char(beamType)
                    case 'paraxial'
                        config.beam.approx_code = 'C';
                    case 'geometric_hat'
                        config.beam.approx_code = 'G';
                    case 'geometric_gaussian'
                        config.beam.approx_code = 'B';
                    otherwise
                        error('Unsupported beam type: %s', beamType);
                end

                % 同步 Bellhop 风格字符串，避免 validate_config_2d 又把 family 覆盖回去
                if ~isfield(config, 'bellhop') || ~isfield(config.bellhop, 'options3') || isempty(config.bellhop.options3)
                    config.bellhop.options3 = 'CCRRR';
                end
                opt3 = upper(char(config.bellhop.options3));
                if numel(opt3) < 5
                    opt3 = [opt3 repmat('R', 1, 5-numel(opt3))];
                    opt3(1) = 'C';  % coherent tl
                end
                opt3(2) = config.beam.approx_code;
                config.bellhop.options3 = opt3;
            case 'step'      % 步长
                config.ray.ds = paramValue;
            case {'integrator', 'integrator_method', 'numerics.integrator'}
                config.numerics.integrator = char(string(paramValue));
            case {'tabular_mode', 'env.profile.tabular_mode'}
                config.env.profile.tabular_mode = char(string(paramValue));
            case {'ssp_interface_jump', 'numerics.ssp_interface_jump'}
                config.numerics.ssp_interface_jump = logical(paramValue);
            case {'ssp_event_tol', 'numerics.ssp_event_tol'}
                config.numerics.ssp_event_tol = double(paramValue);
            case 'cb'        % 海底声速
                config.boundary.bottom.material.c = paramValue;
            case 'rhob'        % 海底密度
                config.boundary.bottom.material.rho = paramValue;
            case 'attnb'        % 海底衰减系数
                config.boundary.bottom.material.attn = paramValue;
            case {'bottom_interp', 'boundary.bottom.interp'}
                config.boundary.bottom.interp = char(string(paramValue));

            case {'bottom_curvature', 'boundary.bottom.curvature_estimator'}
                config.boundary.bottom.curvature_estimator = char(string(paramValue));

            case {'bottom_curvature_on_facet', 'boundary.bottom.use_curvature_on_facet'}
                config.boundary.bottom.use_curvature_on_facet = logical(paramValue);

                % ----------------- 完整路径映射 -----------------
            otherwise
                % 利用 eval 动态解析嵌套的结构体字段名
                % 例如输入 'ray.nbeams'，将执行 config.ray.nbeams = paramValue
                try
                    % 为了防止 eval 解析字符串报错，这里采用安全的结构体赋值策略
                    fieldParts = strsplit(paramName, '.');
                    config = set_nested_field(config, fieldParts, paramValue);
                catch
                    warning('Invalid or unsupported config parameter name: %s', paramName);
                end
        end
    end
end
end

%% 辅助函数：安全地设置嵌套结构体的值
function S = set_nested_field(S, fields, val)
if length(fields) == 1
    S.(fields{1}) = val;
else
    % 如果中间的层级不存在，则初始化为空结构体
    if ~isfield(S, fields{1})
        S.(fields{1}) = struct();
    end
    S.(fields{1}) = set_nested_field(S.(fields{1}), fields(2:end), val);
end
end
