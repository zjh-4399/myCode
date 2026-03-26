function yNext = step_beam_bellhop2(y, h, config)
%STEP_BEAM_BELLHOP2 Bellhop-like 2nd-order beam step.
% 这里采用 midpoint / modified polygon 风格的二阶显式步进：
%   k1   = f(y_n)
%   ymid = y_n + 0.5*h*k1
%   k2   = f(ymid)
%   y_{n+1} = y_n + h*k2
%
% 这不是逐行照抄 Bellhop Fortran 的 Step2D，
% 但在“二阶、两次 RHS、比 RK4 更轻”的意义上是对齐的。

k1 = utils.beam_rhs(y, config);
yMid = y + 0.5 * h * k1;
k2 = utils.beam_rhs(yMid, config);

yNext = y + h * k2;
end