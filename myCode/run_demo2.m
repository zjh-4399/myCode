clc;clear all;
% 声源频率 Hz
f0 = 50;
% 声速剖面的深度向量zw
zw = [0 ;200;250;400;600;800;1000;1200;1400;1600;1800;2000;2200;2400;2600;2800;3000;3200;3400;3600;3800;4000;4200;4400;4600;4800;5000];
% 声速剖面的声速向量cw
cw = [1548.52 ;1530.29;1526.69;1517.78;1509.49;1504.30;1501.38;1500.14;1500.12;1501.02;1502.57;1504.62;1507.02;1509.69;1512.55;1515.56;1518.67;1521.85;1525.10;1528.38;1531.7;1535.04;1538.39;1541.76;1545.14;1548.52;1551.91];
% 声源深度 m
sd = 1400;
% 海底密度 g/cm^3
rhob = 1.8;
% 海底声速 m/s
cb = 1600;
% 海底衰减系数 单位需根据后面的ssp_option(3)设定，目前ssp_option(3)='W'，单位为dB/波长
attnb = 0.8;
% 最远接收距离 km
rmax = 100*1e3;
% 距离方向接收水听器位置设置
fr = 100:100:rmax;
% 深度方向接收水听器位置设置
fd = 5:5:5000;
% 声线出射范围，这里设置的是-90~90
theta=[-90 90];
nbeams = 181;
% % 声线追踪区域深度方向，通常需要比最深接收水听器深一些
% zbox = 5500;%ceil(zmax*1.1);
% % 声线追踪区域距离方向，通常需要比最远接收水听器远一些
% rbox = 101;%ceil(rmax*1.1);
% geometric_gaussian
% geometric_hat
% paraxial
% rk4_fixed  % bellhop2_fixed
config = generate_config('cw', cw, 'zw', zw,...
    'freq',f0,'sd',sd,...
    'fd',fd,'fr',fr,'rr',fr,'rd',1400,...
    'theta',theta,'beam','geometric_gaussian',...
    'cb',cb,'rhob',rhob,'attnb',attnb,...
    'integrator', 'bellhop2_fixed', ...
    'tabular_mode', 'piecewise_linear', ...
    'ssp_interface_jump', true, ...
    'bottom_interp', 'normal_linear', ...
    'bottom_curvature', 'bellhop_dss', ...
    'bottom_curvature_on_facet', false,...
    'nbeams',nbeams, ...
    'output.save_ray', false, ...
    'output.plot_ray', false, ...
    'output.plot_tl', false, ...
    'output.store_beams', false);%, ...
%     'run_type','Ray');
tic
res = run_solver_2d(config);
results = res;
toc
save('./results/results_geometric_gaussian.mat','config','results')

%%
figure
% subplot 132
pcolor(results.field_ranges./1e3, results.field_depths, results.field_tl)
shading flat
axis ij
colormap(flipud(jet))
caxis([60 140])
xlabel('Range(km)')
ylabel('Depth(m)')
title('geometric gaussian')


load('E:\zju_work11\Result_sofar_M1\50Hz\MunkS_ref.shd.mat','pressure')
tl_ref = -20*log10(abs(squeeze(pressure)));
tl_ref(1,:)=[];
tl_ref(:,1)=[];
figure
% subplot 132
pcolor(0.1:0.1:100, 1:1:5e3, tl_ref)
shading flat
axis ij
colormap(flipud(jet))
caxis([60 140])
xlabel('Range(km)')
ylabel('Depth(m)')
title('Scooter')



load('E:\zju_work11\Result_sofar_M1\50Hz\MunkB_CB_M1.shd.mat','pressure')
tl_bellhop = -20*log10(abs(squeeze(pressure)));
tl_bellhop(1,:)=[];
tl_bellhop(:,1)=[];

figure
% subplot 132
pcolor(0.1:0.1:100, 1:1:5e3, tl_bellhop)
shading flat
axis ij
colormap(flipud(jet))
caxis([60 140])
xlabel('Range(km)')
ylabel('Depth(m)')
title('Bellhop')

figure
plot(results.field_ranges./1e3,results.field_tl(1400/5,:),'linewidth',2)
hold on
axis ij
plot(0.1:0.1:100,tl_ref(1400,:),'linewidth',2)
plot(0.1:0.1:100,tl_bellhop(1400,:),'linewidth',2)

xlabel('Range(km)')
ylabel('TL(dB)')
grid on
legend('gbt','WI','Bellhop')
ylim([40 110])
%%
% plot_ray_paths(config, results.beams)
% %%
% figure
% set(gcf,'position',[300 300 1400 300])
% subplot 131
% pcolor(results.geometric_hat.field_ranges./1e3, results.geometric_hat.field_depths, results.geometric_hat.field_tl)
% shading flat
% axis ij
% colormap(flipud(jet))
% caxis([60 140])
% xlabel('Range(km)')
% ylabel('Depth(m)')
% title('geometric hat')
%
% % figure
% subplot 132
% pcolor(results.geometric_gaussian.field_ranges./1e3, results.geometric_gaussian.field_depths, results.geometric_gaussian.field_tl)
% shading flat
% axis ij
% colormap(flipud(jet))
% caxis([60 140])
% xlabel('Range(km)')
% ylabel('Depth(m)')
% title('geometric gaussian')
%
% % figure
% subplot 133
% pcolor(results.paraxial.field_ranges./1e3, results.paraxial.field_depths, results.paraxial.field_tl)
% shading flat
% axis ij
% colormap(flipud(jet))
% caxis([60 140])
% xlabel('Range(km)')
% ylabel('Depth(m)')
% title('paraxial')
