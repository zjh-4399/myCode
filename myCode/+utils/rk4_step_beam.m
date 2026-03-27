function yNext = rk4_step_beam(y, h, config)
%RK4_STEP_BEAM Classical RK4 step for the beam ODE system.
k1 = utils.beam_rhs(y, config);
k2 = utils.beam_rhs(y + 0.5*h*k1, config);
k3 = utils.beam_rhs(y + 0.5*h*k2, config);
k4 = utils.beam_rhs(y + h*k3, config);
yNext = y + (h/6.0) * (k1 + 2*k2 + 2*k3 + k4);
end
