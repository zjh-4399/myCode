function yNext = advance_beam_step(y, h, config)
%ADVANCE_BEAM_STEP Dispatch beam integrator by config.
%
% Supported:
%   rk4_fixed
%   bellhop2_fixed
%
% 优先读取：
%   config.integrator.method
% 若没有，再读取：
%   config.numerics.integrator

% 优化：缓存积分器方法名，消除每步 isfield 检查与字符串转换开销
% 切换 config 时可调用 clear utils.advance_beam_step 重置缓存
persistent cachedMethod
if isempty(cachedMethod)
    cachedMethod = local_get_integrator_name(config);
end

switch cachedMethod
    case {'rk4', 'rk4_fixed'}
        yNext = utils.rk4_step_beam(y, h, config);

    case {'bellhop2', 'bellhop2_fixed', 'midpoint2', 'polygon2'}
        yNext = utils.step_beam_bellhop2(y, h, config);

    otherwise
        error('Unsupported integrator: %s', cachedMethod);
end
end

function method = local_get_integrator_name(config)
method = '';

if isfield(config, 'integrator') && isfield(config.integrator, 'method') ...
        && ~isempty(config.integrator.method)
    method = char(string(config.integrator.method));
elseif isfield(config, 'numerics') && isfield(config.numerics, 'integrator') ...
        && ~isempty(config.numerics.integrator)
    method = char(string(config.numerics.integrator));
else
    method = 'rk4_fixed';
end

method = lower(strtrim(method));
end