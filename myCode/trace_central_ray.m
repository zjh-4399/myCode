function ray = trace_central_ray(config, alpha_deg)
%TRACE_CENTRAL_RAY Trace only the central ray while reusing the full beam propagator.
beam = trace_beam_auxiliary(config, alpha_deg);
ray = struct();
ray.alpha_deg = beam.alpha_deg;
ray.alpha_rad = beam.alpha_rad;
ray.r = beam.r;
ray.z = beam.z;
ray.pr = beam.pr;
ray.pz = beam.pz;
ray.tau = beam.tau;
ray.c = beam.c;
ray.Nr = beam.Nr;
ray.Nz = beam.Nz;
ray.surface_bounces = beam.surface_bounces;
end
