function [qEff, W] = local_geometric_width(config, qReal, c0local, useStent)
dalpha = abs(config.beam.delta_alpha_rad);
if dalpha <= 0
    dalpha = deg2rad(1.0);
end

if nargin < 3 || isempty(c0local) || ~isfinite(c0local)
    c0local = 1500;
end

qAbs = abs(qReal);
Wgeom = config.beam.model.width_scale * qAbs * dalpha / max(c0local, eps);

if useStent && config.beam.model.use_stent
    W = max(Wgeom, config.beam.model.width_floor);
    qFloor = max(c0local * W / max(config.beam.model.width_scale * dalpha, eps), eps);
    if qReal < 0
        qEff = -max(qAbs, qFloor);
    else
        qEff =  max(qAbs, qFloor);
    end
else
    W = Wgeom;
    if qReal == 0
        qEff = eps;
    else
        qEff = qReal;
    end
end
end