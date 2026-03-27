function yOut = reflect_beam_bottom(yIn, hit, config)
%REFLECT_BEAM_BOTTOM Bottom reflection for 2D beam tracing.
% State:
%   y = [r; z; pr; pz; Re(P); Im(P); Re(Q); Im(Q); tau; Re(G); Im(G)]

% 优化：缓存曲率修正系数及平面波反射标志，消除每次反射的多级字段访问和 isfield/strcmpi 检查
% 切换 config 时可调用 clear reflect_beam_bottom 重置缓存
persistent cachedCurvFactor cachedPlaneWave
if isempty(cachedCurvFactor)
    cachedCurvFactor = config.beam.model.curvature_factor;
    cachedPlaneWave  = isfield(config.boundary.bottom, 'reflection') && ...
        strcmpi(config.boundary.bottom.reflection, 'plane_wave');
end

r   = hit.rHit;
z   = hit.zHit;

pr  = yIn(3);
pz  = yIn(4);
P   = complex(yIn(5), yIn(6));
Q   = complex(yIn(7), yIn(8));
tau = yIn(9);
G   = complex(yIn(10), yIn(11));

[cw, cr, cz, ~, ~, ~] = sound_speed_derivatives(config, r, z);

tray = [pr; pz];
nB = hit.nBdry(:) / max(norm(hit.nBdry), eps);

Th = dot(tray, nB);
if Th < 0
    nB = -nB;
    Th = -Th;
end
Th = local_safe_denom(Th);

tB = hit.tBdry(:);
if isempty(tB) || norm(tB) < 1.0e-12
    tB = [nB(2); -nB(1)];
end
tB = tB / max(norm(tB), eps);
Tg = dot(tray, tB);

trayRef = tray - 2.0 * Th * nB;
prNew = trayRef(1);
pzNew = trayRef(2);

gradc = [cr; cz];
rayn = [-tray(2); tray(1)];
cnjump = 2.0 * dot(gradc, rayn);
csjump = 2.0 * dot(gradc, tray);

kappa = 0.0;
if isfield(hit, 'kappa') && ~isempty(hit.kappa) && isfinite(hit.kappa)
    kappa = hit.kappa;
end

RM = Tg / Th;
RN = 2.0 * kappa / (cw^2 * Th);
RN = RN + RM * (2.0 * cnjump - RM * csjump) / cw;
RN = cachedCurvFactor * RN;

PNew = P + Q * RN;
QNew = Q;

Rbot = 1.0 + 0.0i;
if cachedPlaneWave
    hitLocal = hit;
    hitLocal.nBdry = nB;
    hitLocal.tBdry = tB;
    yLocal = yIn;
    yLocal(3) = prNew;
    yLocal(4) = pzNew;
    Rbot = bottom_reflection_coeff_fluid(yLocal, hitLocal, config, cw);
end

GNew = G * Rbot;

yOut = [r; z; prNew; pzNew; ...
        real(PNew); imag(PNew); real(QNew); imag(QNew); ...
        tau; real(GNew); imag(GNew)];
end

function x = local_safe_denom(x)
if abs(x) < 1.0e-10
    x = sign(x + (x == 0)) * 1.0e-10;
end
end
