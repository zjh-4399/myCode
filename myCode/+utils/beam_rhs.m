function dyds = beam_rhs(y, config)
%BEAM_RHS Right-hand side for ray + beam ODE system.
% y = [r; z; pr; pz; Re(P); Im(P); Re(Q); Im(Q); tau; Re(G); Im(G)]

r  = y(1);
z  = y(2);
pr = y(3);
pz = y(4);
PRe = y(5);
PIm = y(6);
QRe = y(7);
QIm = y(8);

[c, cr, cz, crr, crz, czz] = sound_speed_derivatives(config, r, z);

Nr = c * pz;
Nz = -c * pr;
Cnn = crr * Nr^2 + 2.0 * crz * Nr * Nz + czz * Nz^2;

drds   = c * pr;
dzds   = c * pz;
dprds  = -cr / c^2;
dpzds  = -cz / c^2;

dPfac  = -(Cnn / c^2);
dQfac  = c;
dPRe   = dPfac * QRe;
dPIm   = dPfac * QIm;
dQRe   = dQfac * PRe;
dQIm   = dQfac * PIm;

dtauds = 1.0 / c;

dyds = [drds; dzds; dprds; dpzds; ...
        dPRe; dPIm; dQRe; dQIm; ...
        dtauds; 0.0; 0.0];
end
