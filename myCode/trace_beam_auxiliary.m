function [beam, onlineTargets] = trace_beam_auxiliary(config, alpha_deg, onlineTargets, traceOpts)
%TRACE_BEAM_AUXILIARY Trace one beam in 2D.
%
% Optional 3rd input onlineTargets enables segment-level online accumulation:
% every newly accepted ray segment [point_k, point_{k+1}] contributes
% immediately to pressure, instead of waiting for a post-trace sweep over the
% whole beam.
%
% State:
%   y = [r; z; pr; pz; Re(P); Im(P); Re(Q); Im(Q); tau; Re(G); Im(G)]

if nargin < 3
    onlineTargets = [];
end
if nargin < 4 || isempty(traceOpts)
    traceOpts = struct();
end
storeTrace = true;
if isfield(traceOpts, 'storeTrace') && ~isempty(traceOpts.storeTrace)
    storeTrace = logical(traceOpts.storeTrace);
end

alpha_rad = deg2rad(alpha_deg);
[c0, ~, ~, ~, ~, ~] = sound_speed_derivatives(config, config.source.range, config.source.depth);
[P0, Q0, initMeta] = local_initial_conditions(config, c0);

pr0 = cos(alpha_rad) / c0;
pz0 = sin(alpha_rad) / c0;

y = [config.source.range; config.source.depth; pr0; pz0; ...
     real(P0); imag(P0); real(Q0); imag(Q0); 0.0; 1.0; 0.0];

point = local_point_from_state(y, config, c0);
lastPoint = point;
if storeTrace
    nStore = config.integrator.nSteps + 1 + ...
             2 * (config.boundary.surface.max_reflections + config.boundary.bottom.max_reflections) + 1000;
    R   = nan(nStore,1);
    Z   = nan(nStore,1);
    Pr  = nan(nStore,1);
    Pz  = nan(nStore,1);
    Pb  = nan(nStore,1);
    Qb  = nan(nStore,1);
    Tau = nan(nStore,1);
    Gb  = nan(nStore,1);
    idx = 1;
    R(1)   = point.r;
    Z(1)   = point.z;
    Pr(1)  = point.pr;
    Pz(1)  = point.pz;
    Pb(1)  = point.P;
    Qb(1)  = point.Q;
    Tau(1) = point.tau;
    Gb(1)  = point.gain;
else
    R = [];
    Z = [];
    Pr = [];
    Pz = [];
    Pb = [];
    Qb = [];
    Tau = [];
    Gb = [];
    idx = 0;
end

nSurface = 0;
nBottom  = 0;

Aconst = local_source_weight(config, c0, Q0, alpha_rad);
beamInfo = struct( ...
    'Aconst', Aconst, ...
    'alpha_rad', alpha_rad, ...
    'c0', c0, ...
    'cosAlphaSqrt', sqrt(max(cos(alpha_rad), 1.0e-12)), ...
    'gaussianAmpScale', Aconst * (1.0 / sqrt(2.0 * pi)) * sqrt(max(cos(alpha_rad), 1.0e-12)), ...
    'hatAmpScale', Aconst * sqrt(max(cos(alpha_rad), 1.0e-12)), ...
    'q0Source', c0 / max(abs(config.beam.delta_alpha_rad), eps));

for k = 1:config.integrator.nSteps
    ds_left = config.integrator.ds;

    while ds_left > 0
        [eventType, dsHit, botHit, sspHit] = local_predict_event(y, ds_left, config);

        if strcmp(eventType, 'none')
            y = utils.advance_beam_step(y, ds_left, config);
            ds_left = 0.0;
            continue;
        end

        dsHit = min(max(dsHit, 0.0), ds_left);

        if dsHit > 0
            y = utils.advance_beam_step(y, dsHit, config);
            point = local_point_from_state(y, config);
            if ~isempty(onlineTargets)
                onlineTargets = accumulate_segment_to_targets_2d(onlineTargets, lastPoint, point, beamInfo, config);
            end
            lastPoint = point;

            if storeTrace
                idx = idx + 1;
                if idx > numel(R)
                    [R, Z, Pr, Pz, Pb, Qb, Tau, Gb] = local_grow_storage(R, Z, Pr, Pz, Pb, Qb, Tau, Gb);
                end
                R(idx)   = point.r;
                Z(idx)   = point.z;
                Pr(idx)  = point.pr;
                Pz(idx)  = point.pz;
                Pb(idx)  = point.P;
                Qb(idx)  = point.Q;
                Tau(idx) = point.tau;
                Gb(idx)  = point.gain;
            end
        end

        switch eventType
            case 'surface'
                y(2) = 0.0;
                y = reflect_beam_surface(y, config);
                nSurface = nSurface + 1;

            case 'bottom'
                y(1) = botHit.rHit;
                y(2) = botHit.zHit;
                y = reflect_beam_bottom(y, botHit, config);
                nBottom = nBottom + 1;

            case 'ssp'
                y(2) = sspHit.zHit;
                y = utils.apply_ssp_interface_jump_2d(y, sspHit, config);

            otherwise
                error('Unknown eventType: %s', eventType);
        end

        point = local_point_from_state(y, config, lastPoint.c);
        lastPoint = point;

        if storeTrace
            idx = idx + 1;
            if idx > numel(R)
                [R, Z, Pr, Pz, Pb, Qb, Tau, Gb] = local_grow_storage(R, Z, Pr, Pz, Pb, Qb, Tau, Gb);
            end
            R(idx)   = point.r;
            Z(idx)   = point.z;
            Pr(idx)  = point.pr;
            Pz(idx)  = point.pz;
            Pb(idx)  = point.P;
            Qb(idx)  = point.Q;
            Tau(idx) = point.tau;
            Gb(idx)  = point.gain;
        end

        ds_left = ds_left - dsHit;
        if ds_left < 1e-12 * max(1.0, config.integrator.ds)
            ds_left = 0.0;
        end

        if nSurface > config.boundary.surface.max_reflections
            ds_left = 0.0;
        end
        if nBottom > config.boundary.bottom.max_reflections
            ds_left = 0.0;
        end
    end

    point = local_point_from_state(y, config);
    if point.r ~= lastPoint.r || point.z ~= lastPoint.z
        if ~isempty(onlineTargets)
            onlineTargets = accumulate_segment_to_targets_2d(onlineTargets, lastPoint, point, beamInfo, config);
        end
        lastPoint = point;
        if storeTrace
            idx = idx + 1;
            if idx > numel(R)
                [R, Z, Pr, Pz, Pb, Qb, Tau, Gb] = local_grow_storage(R, Z, Pr, Pz, Pb, Qb, Tau, Gb);
            end
            R(idx)   = point.r;
            Z(idx)   = point.z;
            Pr(idx)  = point.pr;
            Pz(idx)  = point.pz;
            Pb(idx)  = point.P;
            Qb(idx)  = point.Q;
            Tau(idx) = point.tau;
            Gb(idx)  = point.gain;
        end
    end

    if ~all(isfinite(y))
        break;
    end

    if abs(point.gain) < config.stop.min_gain
        break;
    end

    if y(1) > config.stop.max_range || y(2) > config.stop.max_depth
        break;
    end
end

if storeTrace
    R   = R(1:idx);
    Z   = Z(1:idx);
    Pr  = Pr(1:idx);
    Pz  = Pz(1:idx);
    Pb  = Pb(1:idx);
    Qb  = Qb(1:idx);
    Tau = Tau(1:idx);
    Gb  = Gb(1:idx);

    [c, ~, ~, crr, crz, czz] = sound_speed_derivatives(config, R, Z);
    Nr = c .* Pz;
    Nz = -c .* Pr;
    Cnn = crr .* Nr.^2 + 2.0 * crz .* Nr .* Nz + czz .* Nz.^2;
else
    c = [];
    Nr = [];
    Nz = [];
    Cnn = [];
end

beam = struct();
beam.alpha_deg = alpha_deg;
beam.alpha_rad = alpha_rad;
beam.r = R;
beam.z = Z;
beam.pr = Pr;
beam.pz = Pz;
beam.P = Pb;
beam.Q = Qb;
beam.tau = Tau;
beam.c = c;
beam.Nr = Nr;
beam.Nz = Nz;
beam.Cnn = Cnn;
beam.P0 = P0;
beam.Q0 = Q0;
beam.init = initMeta;
beam.Aconst = Aconst;
beam.surface_bounces = nSurface;
beam.gain = Gb;
beam.bottom_bounces = nBottom;
end


function point = local_point_from_state(y, config, c)
if nargin < 3 || isempty(c)
    [c, ~, ~, ~, ~, ~] = sound_speed_derivatives(config, y(1), y(2));
end
point = struct();
point.r = y(1);
point.z = y(2);
point.pr = y(3);
point.pz = y(4);
point.P = complex(y(5), y(6));
point.Q = complex(y(7), y(8));
point.tau = y(9);
point.gain = complex(y(10), y(11));
point.c = c;
point.Nr = c * y(4);
point.Nz = -c * y(3);
end

% =========================================================================
% Local helpers
% =========================================================================

function [eventType, dsHit, botHit, sspHit] = local_predict_event(y, dsMax, config)
eventType = 'none';
dsHit = inf;
botHit = local_empty_bottom_hit();
sspHit = local_empty_ssp_hit();

[c, ~, ~, ~, ~, ~] = sound_speed_derivatives(config, y(1), y(2));
drds = c * y(3);
dzds = c * y(4);

yPred = y;
yPred(1) = y(1) + drds * dsMax;
yPred(2) = y(2) + dzds * dsMax;

zTol = 1e-12 * max(1.0, abs(y(2)));
sTol = 1.0e-14 * max(1.0, dsMax);

if isfield(config, 'boundary') && isfield(config.boundary, 'surface') ...
        && config.boundary.surface.enabled

    if (y(2) > zTol) && (dzds < -sTol)
        dsSurf = -y(2) / dzds;
        if (dsSurf > sTol) && (dsSurf <= dsMax + sTol)
            eventType = 'surface';
            dsHit = dsSurf;
        end
    end
end

botCand = local_find_bottom_hit(y, yPred, dsMax, config);
if botCand.exists && botCand.dsHit < dsHit
    eventType = 'bottom';
    dsHit = botCand.dsHit;
    botHit = botCand;
end

sspCand = local_find_ssp_interface_hit(y, dsMax, config, c);
if sspCand.exists && sspCand.dsHit < dsHit
    eventType = 'ssp';
    dsHit = sspCand.dsHit;
    sspHit = sspCand;
end

if ~isfinite(dsHit)
    eventType = 'none';
    dsHit = inf;
end
end

function hit = local_find_ssp_interface_hit(y, dsMax, config, c0)
hit = local_empty_ssp_hit();

if ~local_ssp_jump_enabled(config)
    return;
end
if ~strcmpi(config.env.profile.type, 'tabular')
    return;
end
if ~strcmpi(local_get_tabular_mode(config), 'piecewise_linear')
    return;
end

zw = config.env.profile.zw(:);
if numel(zw) < 3
    return;
end

if nargin < 4 || isempty(c0)
    [c0, ~, ~, ~, ~, ~] = sound_speed_derivatives(config, y(1), y(2));
end
dzds = c0 * y(4);
if abs(dzds) < 1.0e-14 * max(1.0, dsMax)
    return;
end

z0 = y(2);
internalZ = zw(2:end-1);

zTol = local_ssp_event_tol(config) * max(1.0, max(abs([z0; internalZ])));
sTol = 1.0e-14 * max(1.0, dsMax);

if dzds > 0
    cand = internalZ(internalZ > z0 + zTol);
    if isempty(cand)
        return;
    end
    zIf = cand(1);
else
    cand = internalZ(internalZ < z0 - zTol);
    if isempty(cand)
        return;
    end
    zIf = cand(end);
end

dsIf = (zIf - z0) / dzds;
if ~(dsIf > sTol && dsIf <= dsMax + sTol)
    return;
end

j = find(abs(zw - zIf) <= zTol, 1, 'first');
if isempty(j)
    return;
end

hit.exists = true;
hit.dsHit = dsIf;
hit.zHit = zIf;
hit.interfaceIndex = j;
end

function hit = local_empty_bottom_hit()
hit = struct( ...
    'exists', false, ...
    'rHit',   nan, ...
    'zHit',   nan, ...
    'segId',  nan, ...
    'dsHit',  inf, ...
    'nBdry',  [], ...
    'tBdry',  [], ...
    'kappa',  [] );
end

function hit = local_empty_ssp_hit()
hit = struct( ...
    'exists', false, ...
    'dsHit', inf, ...
    'zHit', nan, ...
    'interfaceIndex', nan );
end

function tf = local_ssp_jump_enabled(config)
tf = true;
if isfield(config, 'numerics') && isfield(config.numerics, 'ssp_interface_jump') ...
        && ~isempty(config.numerics.ssp_interface_jump)
    tf = logical(config.numerics.ssp_interface_jump);
end
end

function tol = local_ssp_event_tol(config)
tol = 1.0e-10;
if isfield(config, 'numerics') && isfield(config.numerics, 'ssp_event_tol') ...
        && ~isempty(config.numerics.ssp_event_tol)
    tol = double(config.numerics.ssp_event_tol);
end
end

function mode = local_get_tabular_mode(config)
mode = 'piecewise_linear';
if isfield(config, 'env') && isfield(config.env, 'profile') ...
        && isfield(config.env.profile, 'tabular_mode') ...
        && ~isempty(config.env.profile.tabular_mode)
    mode = lower(strtrim(char(string(config.env.profile.tabular_mode))));
end
end

function [R, Z, Pr, Pz, Pb, Qb, Tau, Gb] = local_grow_storage(R, Z, Pr, Pz, Pb, Qb, Tau, Gb)
nOld = numel(R);
nNew = nOld + max(1000, ceil(0.5 * nOld));

R(nNew)   = nan;
Z(nNew)   = nan;
Pr(nNew)  = nan;
Pz(nNew)  = nan;
Pb(nNew)  = nan;
Qb(nNew)  = nan;
Tau(nNew) = nan;
Gb(nNew)  = nan;
end

function [P0, Q0, meta] = local_initial_conditions(config, c0)
beamModel = config.beam.model;
dalpha = abs(config.beam.delta_alpha_rad);
if dalpha <= 0
    dalpha = deg2rad(1.0);
end

epsSpace = 2.0 * c0^2 / (config.omega * dalpha^2);
epsSpace = epsSpace * beamModel.epmult;

switch beamModel.family
    case 'paraxial'
        switch beamModel.init_mode
            case {'cerveny_space_filling','space_filling'}
                epsilon0 = epsSpace;
            case 'fixed_width'
                if isempty(beamModel.fixed_width)
                    error('beam.init_code = F requires config.beam.fixed_width for paraxial beams.');
                end
                epsilon0 = 0.5 * config.omega * beamModel.fixed_width.^2;
            case 'minimum_width'
                epsilonFloor = 0.5 * config.omega * beamModel.width_floor.^2;
                epsilon0 = max(epsSpace, epsilonFloor);
            case 'wkb_like'
                epsilon0 = 0.5 * epsSpace;
            otherwise
                error('Unsupported paraxial beam init mode: %s', beamModel.init_mode);
        end

        if ~isempty(beamModel.manual_epsilon)
            epsilon0 = beamModel.manual_epsilon;
        end

        P0 = 1.0 + 0i;
        Q0 = 1i * epsilon0;

        meta = struct();
        meta.mode = beamModel.init_mode;
        meta.epsilon0 = epsilon0;
        meta.epsilon_space = epsSpace;

    case {'geometric_hat','geometric_gaussian'}
        P0 = 1.0 + 0i;
        Q0 = 0.0 + 0i;
        meta = struct();
        meta.mode = 'real_ray_tube';
        meta.epsilon0 = 0.0;
        meta.epsilon_space = epsSpace;

    otherwise
        error('Unsupported beam family in initial condition builder: %s', beamModel.family);
end
end

function Aconst = local_source_weight(config, c0, Q0, alpha_rad)
sourcePattern = 1.0;
if isfield(config.source, 'beam_pattern') && strcmpi(config.source.beam_pattern, 'omni')
    sourcePattern = 1.0;
end

switch config.beam.model.family
    case 'paraxial'
        Aconst = sourcePattern * (1.0 / c0) * exp(1i*pi/4) * ...
            sqrt(Q0 * config.omega * max(cos(alpha_rad), 1e-12) / (2.0*pi));
    case {'geometric_hat','geometric_gaussian'}
        Aconst = sourcePattern;
    otherwise
        error('Unsupported beam family: %s', config.beam.model.family);
end
end
