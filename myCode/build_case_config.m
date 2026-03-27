function config = build_case_config(case_id)
%BUILD_CASE_CONFIG Compatibility wrapper around CREATE_BELLHOP_STYLE_CONFIG.
% This keeps the historical file name while moving user-facing parameters
% into a Bellhop-style config tree.

if nargin < 1 || isempty(case_id)
    case_id = 'munk';
end
config = create_bellhop_style_config(case_id);
config = validate_config_2d(config);
end
