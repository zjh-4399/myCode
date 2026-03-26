function beamModel = parse_beam_options_2d(config)
%PARSE_BEAM_OPTIONS_2D Parse Bellhop-style beam options into solver options.
%
% Supported family / approx-code combinations in this 2D project:
%   paraxial           <-> C
%   geometric_hat      <-> G
%   geometric_gaussian <-> B
%
% Bellhop's ray-centered option (R) is parsed but not implemented.

beam = config.beam;

approx_code = upper(local_char_or_default(beam, 'approx_code', 'AUTO'));
family = lower(local_char_or_default(beam, 'family', 'auto'));
init_code = upper(local_char_or_default(beam, 'init_code', 'C'));
curv_code = upper(local_char_or_default(beam, 'curvature_code', 'S'));

if strcmpi(family, 'auto')
    switch approx_code
        case 'C'
            family = 'paraxial';
        case 'G'
            family = 'geometric_hat';
        case 'B'
            family = 'geometric_gaussian';
        case {'AUTO','A'}
            family = 'paraxial';
            approx_code = 'C';
        case 'R'
            error(['approx_code = ''R'' is reserved for ray-centered beam influence, but it is ', ...
                'not implemented in this 2D MATLAB project.']);
        otherwise
            error('Unsupported beam.approx_code: %s', approx_code);
    end
else
    switch family
        case 'paraxial'
            if any(strcmp(approx_code, {'AUTO','A'}))
                approx_code = 'C';
            elseif ~strcmp(approx_code, 'C')
                error('family = paraxial only supports approx_code = C in this project.');
            end
        case 'geometric_hat'
            if any(strcmp(approx_code, {'AUTO','A'}))
                approx_code = 'G';
            elseif ~strcmp(approx_code, 'G')
                error('family = geometric_hat only supports approx_code = G in this project.');
            end
        case 'geometric_gaussian'
            if any(strcmp(approx_code, {'AUTO','A'}))
                approx_code = 'B';
            elseif ~strcmp(approx_code, 'B')
                error('family = geometric_gaussian only supports approx_code = B in this project.');
            end
        otherwise
            error('Unsupported beam.family: %s', family);
    end
end

switch init_code
    case 'C'
        init_mode = 'cerveny_space_filling';
    case 'F'
        init_mode = 'space_filling';
    case 'M'
        init_mode = 'minimum_width';
    case 'W'
        init_mode = 'wkb_like';
    otherwise
        error('Unsupported beam.init_code: %s', init_code);
end

switch curv_code
    case 'S'
        curvature_factor = 1.0;
        curvature_mode = 'standard';
    case 'D'
        curvature_factor = 2.0;
        curvature_mode = 'doubling';
    case 'Z'
        curvature_factor = 0.0;
        curvature_mode = 'zeroing';
    otherwise
        error('Unsupported beam.curvature_code: %s', curv_code);
end

beamModel = struct();
beamModel.family = family;
beamModel.approx_code = approx_code;
beamModel.init_code = init_code;
beamModel.init_mode = init_mode;
beamModel.curvature_code = curv_code;
beamModel.curvature_mode = curvature_mode;
beamModel.curvature_factor = curvature_factor;
beamModel.shift_enabled = logical(beam.shift_enabled);
beamModel.width_floor = beam.width_floor;
beamModel.width_scale = beam.width_scale;
beamModel.window_radii = beam.window_radii;
beamModel.ibwin = beam.ibwin;
beamModel.epmult = beam.epmult;
beamModel.fixed_width = beam.fixed_width;
beamModel.manual_epsilon = beam.manual_epsilon;
beamModel.auto_window = logical(beam.auto_window);
beamModel.use_stent = logical(beam.use_stent);
beamModel.name = sprintf('%s/%s/%s/%s', family, approx_code, init_code, curv_code);
end

function value = local_char_or_default(s, fieldName, defaultValue)
if isfield(s, fieldName) && ~isempty(s.(fieldName))
    value = string(s.(fieldName));
    value = char(value);
else
    value = defaultValue;
end
end
