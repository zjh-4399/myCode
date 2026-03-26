function [pressure, beams] = accumulate_pressure_from_all_beams(config, receiver_depths, receiver_ranges, beams)
%ACCUMULATE_PRESSURE_FROM_ALL_BEAMS Sum complex pressure from all beams.

if nargin < 4 || isempty(beams)
    beams = cell(config.ray.nbeams, 1);
    for ib = 1:config.ray.nbeams
        beams{ib} = trace_beam_auxiliary(config, config.ray.theta_list_deg(ib));
    end
end

nZ = numel(receiver_depths);
nR = numel(receiver_ranges);
pressure = complex(zeros(nZ, nR));

for iz = 1:nZ
    zr = receiver_depths(iz);
    for ib = 1:numel(beams)
        pressure(iz, :) = pressure(iz, :) + ...
            interpolate_beam_to_receiver(beams{ib}, zr, receiver_ranges, config);
    end
end
end


% function config = apply_bellhop_option_strings_2d(config)
% %APPLY_BELLHOP_OPTION_STRINGS_2D Apply Bellhop-like OPTIONS3/OPTIONS4 strings.
% % This project supports only a 2D subset.
% %
% % OPTIONS3 positions used here:
% %   1: output type      R/C/A/E/I/S
% %   2: beam approximation G/C/R/B
% %   3: beam shift       <blank>/S/*
% %   4: source type      R/X
% %   5: receiver grid    R/I
% %
% % OPTIONS4 positions used here:
% %   1: beam init        C/F/M/W
% %   2: curvature        D/S/Z
% %
% % Unsupported combinations raise explicit errors.
% 
% if ~isfield(config, 'bellhop')
%     return;
% end
% 
% if ~isfield(config.bellhop, 'options3') || isempty(config.bellhop.options3)
%     return;
% end
% 
% opt3 = upper(char(string(config.bellhop.options3)));
% opt3 = opt3(:).';
% if numel(opt3) < 5
%     opt3 = [opt3 repmat(' ', 1, 5 - numel(opt3))];
% else
%     opt3 = opt3(1:5);
% end
% 
% if ~isfield(config.bellhop, 'options4') || isempty(config.bellhop.options4)
%     opt4 = '  ';
% else
%     opt4 = upper(char(string(config.bellhop.options4)));
%     opt4 = opt4(:).';
%     if numel(opt4) < 2
%         opt4 = [opt4 repmat(' ', 1, 2 - numel(opt4))];
%     else
%         opt4 = opt4(1:2);
%     end
% end
% 
% % OPTIONS3(1): output type
% switch opt3(1)
%     case {'R','C','A','E','I','S'}
%         config.output.type = opt3(1);
%     case ' '
%         % leave as-is
%     otherwise
%         error('Unsupported Bellhop OPTIONS3(1): %s', opt3(1));
% end
% 
% % OPTIONS3(2): approximation / beam family
% switch opt3(2)
%     case ' '
%         % leave as-is
%     case 'G'
%         config.beam.family = 'geometric_hat';
%         config.beam.approx_code = 'G';
%     case 'C'
%         config.beam.family = 'paraxial';
%         config.beam.approx_code = 'C';
%     case 'B'
%         config.beam.family = 'geometric_gaussian';
%         config.beam.approx_code = 'B';
%     case 'R'
%         error(['Bellhop OPTIONS3(2) = R (ray-centered beam influence) is not implemented ', ...
%             'in this 2D MATLAB project.']);
%     otherwise
%         error('Unsupported Bellhop OPTIONS3(2): %s', opt3(2));
% end
% 
% % OPTIONS3(3): shift / legacy no-shift code
% switch opt3(3)
%     case {' ','R'}
%         config.beam.shift_enabled = false;
%     case 'S'
%         config.beam.shift_enabled = true;
%     case '*'
%         error('Bellhop source beam pattern files (*) are not implemented in this MATLAB project.');
%     otherwise
%         error(['Unsupported Bellhop OPTIONS3(3): %s. ', ...
%             'Use blank or R for no shift, S for shift.'], opt3(3));
% end
% 
% % OPTIONS3(4): source type
% switch opt3(4)
%     case {' ','R'}
%         config.source.coordinate_type = 'cylindrical_point';
%     case 'X'
%         error('Line-source Cartesian source (OPTIONS3(4)=X) is not implemented in this 2D MATLAB project.');
%     otherwise
%         error('Unsupported Bellhop OPTIONS3(4): %s', opt3(4));
% end
% 
% % OPTIONS3(5): receiver grid
% switch opt3(5)
%     case {' ','R'}
%         config.receiver.grid_type = 'rectilinear';
%     case 'I'
%         config.receiver.grid_type = 'irregular';
%     otherwise
%         error('Unsupported Bellhop OPTIONS3(5): %s', opt3(5));
% end
% 
% % OPTIONS4(1): init type
% switch opt4(1)
%     case ' '
%         % leave as-is
%     case {'C','F','M','W'}
%         config.beam.init_code = opt4(1);
%     otherwise
%         error('Unsupported Bellhop OPTIONS4(1): %s', opt4(1));
% end
% 
% % OPTIONS4(2): curvature
% switch opt4(2)
%     case ' '
%         % leave as-is
%     case {'D','S','Z'}
%         config.beam.curvature_code = opt4(2);
%     otherwise
%         error('Unsupported Bellhop OPTIONS4(2): %s', opt4(2));
% end
% 
% if isfield(config.bellhop, 'component') && ~isempty(config.bellhop.component)
%     config.output.component = char(string(config.bellhop.component));
% end
% 
% config.bellhop.options3_applied = opt3;
% config.bellhop.options4_applied = opt4;
% end
% 
% 
