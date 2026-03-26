function config = create_bellhop_style_config(preset_name)
%CREATE_BELLHOP_STYLE_CONFIG Create a 2D beam-tracing config with Bellhop-style options.
%   config = CREATE_BELLHOP_STYLE_CONFIG()
%   config = CREATE_BELLHOP_STYLE_CONFIG('munk')
%
% This file exposes user-facing configuration only. Use VALIDATE_CONFIG_2D
% before solving so all defaults, aliases, and parsers are applied.
%
% Supported beam families in this 2D project:
%   - paraxial
%   - geometric_hat
%   - geometric_gaussian
%
% Supported output types:
%   - ray / R
%   - coherent_tl / C
%
% Not implemented:
%   - arrivals / A
%   - eigenray / E
%   - incoherent_tl / I
%   - semicoherent_tl / S

if nargin < 1 || isempty(preset_name)
    preset_name = 'munk';
end
preset_name = lower(strtrim(preset_name));

% ----------------------------
% Environment
% ----------------------------
config.reference = '2D beam tracing starter with Bellhop-style parameter interface';
config.dimension = '2D';
config.env.profile.type = 'munk';
config.env.profile.params.c0 = 1500.0;
config.env.profile.params.epsilon = 0.00737;
config.env.profile.params.zc = 1300.0;
config.env.profile.params.z_flat = 5000.0;
config.env.profile.tabular_mode = 'piecewise_linear';

% ----------------------------
% Source / receiver
% ----------------------------
config.source.range = 0.0;
config.source.depth = 1000.0;
config.source.freq_hz = 50.0;
config.source.beam_pattern = 'omni';

config.receiver.depths = 800.0;
config.receiver.ranges = linspace(0.25, 100000.0, 800);
config.receiver.grid_type = 'rectilinear';

% ----------------------------
% Ray fan
% ----------------------------
config.ray.nbeams = 88;
config.ray.auto_nbeams = false;
config.ray.theta_min_deg = -14.66;
config.ray.theta_max_deg = 14.66;
config.ray.theta_list_deg = [];
config.ray.ds = 200.0;
config.ray.zmax = 5000.0;
config.ray.rmax = 100000.0;

% ----------------------------
% Output
% ----------------------------
config.output.type = 'coherent_tl';   % 'ray' or 'coherent_tl'
config.output.component = 'pressure';
config.output.save_ray = false;
config.output.save_field = false;
config.output.plot_ray = true;
config.output.plot_tl = true;
config.output.store_beams = false;
config.output.ray_max_points = 500000;
config.output.ray_file_name = 'beam_trace.ray';
config.output.results_dir = fullfile(pwd, 'results');
config.output.fig_dir = fullfile(pwd, 'figures');

% Optional field grid
config.field.enabled = false;
config.field.depths = linspace(0.0, 5000.0, 251);
config.field.ranges = linspace(0.25, 100000.0, 401);

% ----------------------------
% Beam control
% ----------------------------
config.beam.family = 'paraxial';      % paraxial / geometric_hat / geometric_gaussian / auto
config.beam.approx_code = 'C';        % C/G/B or 'auto'
config.beam.init_code = 'C';          % C/F/M/W
config.beam.curvature_code = 'S';     % S/D/Z
config.beam.shift_enabled = false;
config.beam.width_floor = 5.0;
config.beam.width_scale = 1.0;
config.beam.fixed_width = [];
config.beam.manual_epsilon = [];
config.beam.epmult = 1.0;
config.beam.rloop = 1;
config.beam.isingl = 0;
config.beam.nimage = 0;
config.beam.ibwin = 5;
config.beam.window_radii = 5.0;
config.beam.auto_window = true;
config.beam.use_stent = true;

% ----------------------------
% Boundaries
% ----------------------------
config.boundary.surface.type = 'pressure_release';
config.boundary.surface.model = 'grazing_empirical';   % basic / grazing_empirical
config.boundary.surface.enabled = true;
config.boundary.surface.max_reflections = 100;

config.boundary.bottom.type = 'none';
config.boundary.bottom.model = 'none';
config.boundary.bottom.enabled = false;
config.boundary.bottom.interp = 'normal_linear';
config.boundary.bottom.curvature_estimator = 'bellhop_dss';
config.boundary.bottom.use_curvature_on_facet = false;

% ----------------------------
% Numerics
% ----------------------------
config.numerics.integrator = 'rk4_fixed';   % 'rk4_fixed' / 'bellhop2_fixed'
config.numerics.branch_tracking = 'continuous_sqrt';
config.numerics.interp_scheme = 'linear';
config.numerics.max_steps = [];
config.numerics.store_all_steps = true;

config.numerics.ssp_interface_jump = true;
config.numerics.ssp_event_tol = 1e-10;

% Bellhop-style early termination when cumulative reflection amplitude is tiny
config.stop.min_gain = 5.0e-3;

% ----------------------------
% Bellhop-style string interface (optional)
% ----------------------------
config.bellhop.options3 = '';     % e.g. 'CGRRR', 'CBRRR', 'CCRRR'
config.bellhop.options4 = '';     % e.g. 'MS', 'CZ', 'FD'
config.bellhop.component = '';    % only meaningful for ray-centered beams, not implemented here

% ----------------------------
% Plot / metadata
% ----------------------------
config.title = '2D beam tracing case';
config.notes = '';
config.units.length = 'm';
config.units.length_name = 'meter';
config.units.speed = 'm/s';

switch preset_name
    case {'munk','case3','munk_profile'}
        config.title = 'Munk profile';
        config.env.profile.type = 'munk';
        config.env.profile.params.c0 = 1500.0;
        config.env.profile.params.epsilon = 0.00737;
        config.env.profile.params.zc = 1300.0;
        config.env.profile.params.z_flat = 5000.0;

        config.source.depth = 1000.0;
        config.source.freq_hz = 50.0;
        config.receiver.depths = 800.0;
        config.receiver.ranges = linspace(125, 100000.0, 800);
        config.ray.nbeams = 88;
        config.ray.theta_min_deg = -14.66;
        config.ray.theta_max_deg = 14.66;
        config.ray.ds = 500.0;
        config.ray.rmax = 100000.0;
        config.ray.zmax = 5000.0;
        config.output.type = 'coherent_tl';
        config.field.enabled = true;
        config.field.depths = linspace(0.0, 5000.0, 251);
        config.field.ranges = linspace(250, 100000.0, 400);
        config.boundary.surface.model = 'grazing_empirical';
        config.beam.family = 'paraxial';
        config.beam.approx_code = 'C';
        config.bellhop.options3 = 'CCRRR';
        config.bellhop.options4 = 'CS';
        config.notes = 'Preset approximating the published 2D Munk benchmark.';

    case {'downward_shallow','case1','shallow'}
        config.title = 'Downward refracting profile - shallow source';
        config.units.length = 'yd';
        config.units.length_name = 'yard';
        config.units.speed = 'yd/s';

        config.env.profile.type = 'downward';
        config.env.profile.params.c0 = 1677.3319;
        config.env.profile.params.gamma = -1.2286762;

        config.source.depth = 66.7;
        config.source.freq_hz = 2000.0;
        config.receiver.depths = 66.7;
        config.receiver.ranges = linspace(500.0, 1000.0, 501);

        config.ray.nbeams = 101;
        config.ray.theta_min_deg = -25.0;
        config.ray.theta_max_deg = 0.0;
        config.ray.ds = 10.0;
        config.ray.rmax = 1200.0;
        config.ray.zmax = 5000.0;

        config.output.type = 'coherent_tl';
        config.field.enabled = false;
        config.boundary.surface.model = 'basic';
        config.beam.family = 'paraxial';
        config.beam.approx_code = 'C';
        config.bellhop.options3 = 'CCRRR';
        config.bellhop.options4 = 'CS';

    case {'downward_deep','case2','deep'}
        config.title = 'Downward refracting profile - deep source';
        config.units.length = 'yd';
        config.units.length_name = 'yard';
        config.units.speed = 'yd/s';

        config.env.profile.type = 'downward';
        config.env.profile.params.c0 = 1677.3319;
        config.env.profile.params.gamma = -1.2286762;

        config.source.depth = 1000.0;
        config.source.freq_hz = 2000.0;
        config.receiver.depths = 800.0;
        config.receiver.ranges = linspace(3100.0, 3170.0, 281);

        config.ray.nbeams = 401;
        config.ray.theta_min_deg = -60.0;
        config.ray.theta_max_deg = -20.0;
        config.ray.ds = 10.0;
        config.ray.rmax = 3500.0;
        config.ray.zmax = 5000.0;

        config.output.type = 'coherent_tl';
        config.field.enabled = false;
        config.boundary.surface.model = 'basic';
        config.beam.family = 'paraxial';
        config.beam.approx_code = 'C';
        config.bellhop.options3 = 'CCRRR';
        config.bellhop.options4 = 'CS';

    otherwise
        error('Unknown preset_name: %s', preset_name);
end
end
