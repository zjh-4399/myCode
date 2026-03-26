function outModel = parse_output_options_2d(config)
%PARSE_OUTPUT_OPTIONS_2D Parse Bellhop-style output options.
% Supported outputs:
%   R / ray          -> ray-only tracing and plotting
%   C / coherent_tl  -> coherent complex pressure and TL
%
% Explicitly unsupported:
%   A / arrivals
%   E / eigenray
%   I / incoherent_tl
%   S / semicoherent_tl

typeRaw = config.output.type;
if isstring(typeRaw); typeRaw = char(typeRaw); end
typeRaw = lower(strtrim(typeRaw));

outModel = struct();
outModel.plot_ray = logical(config.output.plot_ray);
outModel.plot_tl = logical(config.output.plot_tl);
outModel.save_ray = logical(config.output.save_ray);
outModel.save_field = logical(config.output.save_field);
outModel.component = lower(strtrim(config.output.component));

switch typeRaw
    case {'r','ray','raytrace'}
        outModel.type = 'ray';
        outModel.type_code = 'R';
    case {'c','coherent_tl','coherent','tl'}
        outModel.type = 'coherent_tl';
        outModel.type_code = 'C';
    case {'a','arrivals'}
        error('output.type = arrivals/A is not implemented in this 2D MATLAB project.');
    case {'e','eigenray'}
        error('output.type = eigenray/E is not implemented in this 2D MATLAB project.');
    case {'i','incoherent_tl'}
        error('output.type = incoherent_tl/I is not implemented in this 2D MATLAB project.');
    case {'s','semicoherent_tl'}
        error('output.type = semicoherent_tl/S is not implemented in this 2D MATLAB project.');
    otherwise
        error('Unsupported output.type: %s', typeRaw);
end

if ~strcmp(outModel.component, 'pressure')
    error('Only output.component = pressure is implemented in this 2D MATLAB project.');
end
end
