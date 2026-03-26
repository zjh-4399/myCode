restoredefaultpath;
rehash toolboxcache;
cd(''E:/zju_work11/Munk_M_0'');
addpath(''myCode'');

config = generate_config(''cw'', [1500;1490;1480;1490;1500], ...
    ''zw'', [0;500;1000;1500;2000], ...
    ''freq'', 50, ''sd'', 1000, ...
    ''fd'', 0:20:2000, ''fr'', 500:500:20000, ...
    ''rr'', 500:500:20000, ''rd'', 1000, ...
    ''theta'', [-25 25], ''nbeams'', 41, ...
    ''beam'', ''geometric_gaussian'', ...
    ''integrator'', ''bellhop2_fixed'', ...
    ''tabular_mode'', ''piecewise_linear'', ...
    ''output.save_ray'', false, ''output.plot_ray'', false, ...
    ''output.plot_tl'', false, ''output.store_beams'', false);

profile clear;
profile on -timer real;
res = run_solver_2d(config); %#ok<NASGU>
profile off;

p = profile(''info'');
T = struct2table(p.FunctionTable);
if ~isempty(T)
    T = sortrows(T, ''TotalTime'', ''descend'');
end

outFile = 'E:/zju_work11/Munk_M_0/myCode/results/profile_top.txt';
fid = fopen(outFile, 'w');
if fid < 0
    error('Cannot open profile_top.txt for writing');
end
fprintf(fid, 'nFunctions=%d\n', height(T));
N = min(20, height(T));
for i = 1:N
    fn = T.FunctionName{i};
    tt = T.TotalTime(i);
    nc = T.NumCalls(i);
    fprintf(fid, '%02d | %10.6f s | %8d calls | %s\n', i, tt, nc, fn);
end
fclose(fid);

save('E:/zju_work11/Munk_M_0/myCode/results/profile_workspace.mat', 'p');
