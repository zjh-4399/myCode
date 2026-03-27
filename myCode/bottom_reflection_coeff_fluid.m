function R = bottom_reflection_coeff_fluid(yIn, hit, config, cw)
%BOTTOM_REFLECTION_COEFF_FLUID Plane-wave reflection coefficient for fluid half-space.

% 优化：缓存底部材料参数及预计算复数声速 cbx，消除每次反弹的多级字段访问和重复计算
% 切换 config 时可调用 clear bottom_reflection_coeff_fluid 重置缓存
persistent cachedRho2 cachedCbx
if isempty(cachedRho2)
    cachedRho2 = config.boundary.bottom.material.rho;
    cachedCbx  = complex_bottom_speed_from_dbwl( ...
        config.boundary.bottom.material.c, ...
        config.boundary.bottom.material.attn);
end

pr = yIn(3);
pz = yIn(4);

tB = hit.tBdry(:) / max(norm(hit.tBdry), eps);

% 切向慢度守恒
pPar = abs(pr * tB(1) + pz * tB(2));

rho1 = 1.0;  % 海水，若你后面加入 config.env.rho_water，就替换这里
rho2 = cachedRho2;
cbx  = cachedCbx;

q1 = sqrt(complex(1.0 / cw^2 - pPar^2, 0.0));
q2 = sqrt((1.0 / cbx)^2 - pPar^2);

% 选取向下/衰减的透射分支
if imag(q2) < 0
    q2 = -q2;
end

R = (rho2 * q1 - rho1 * q2) / (rho2 * q1 + rho1 * q2);
end

function cbx = complex_bottom_speed_from_dbwl(cb, attn_db_per_wl)
% attn: dB / wavelength
eta = attn_db_per_wl * log(10.0) / (40.0 * pi);
cbx = cb / (1.0 + 1i * eta);
end