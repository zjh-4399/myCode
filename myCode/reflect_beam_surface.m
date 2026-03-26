function yOut = reflect_beam_surface(yIn, config)
%REFLECT_BEAM_SURFACE Bellhop-consistent pressure-release surface reflection.
% State:
%   y = [r; z; pr; pz; Re(P); Im(P); Re(Q); Im(Q); tau; Re(G); Im(G)]

r   = yIn(1);
pr  = yIn(3);
pz  = yIn(4);
P   = complex(yIn(5), yIn(6));
Q   = complex(yIn(7), yIn(8));
tau = yIn(9);
G   = complex(yIn(10), yIn(11));

[cSurf, crSurf, czSurf, ~, ~, ~] = sound_speed_derivatives(config, r, 0.0);

tray = [pr; pz];
tB = [1.0; 0.0];
nB = [0.0; -1.0];

Tg = dot(tray, tB);
Th = dot(tray, nB);
if abs(Th) < 1.0e-10
    Th = sign(Th + (Th == 0)) * 1.0e-10;
end

gradc = [crSurf; czSurf];
rayn = [-tray(2); tray(1)];
cnjump = 2.0 * dot(gradc, rayn);
csjump = 2.0 * dot(gradc, tray);

% Bellhop flips the normal-gradient term and curvature term for top reflection.
cnjump = -cnjump;
RM = Tg / Th;
RN = RM * (2.0 * cnjump - RM * csjump) / cSurf;
RN = -RN;
RN = config.beam.model.curvature_factor * RN;

prNew = pr;
pzNew = -pz;
PNew = P + Q * RN;
QNew = Q;
GNew = -G;

yOut = [r; 0.0; prNew; pzNew; ...
        real(PNew); imag(PNew); real(QNew); imag(QNew); ...
        tau; real(GNew); imag(GNew)];
end
