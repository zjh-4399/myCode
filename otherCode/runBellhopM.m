figure
pcolor(r_Ray./1e3,z_Ray,TL_Ray)
axis ij
shading flat
title('Ray')
colormap(flipud(jet))
caxis([60 110])
xlabel('距离(km)')
ylabel('深度(m)')
c1 = colorbar;
ylabel(c1,'TL(dB)')


figure
pcolor(r_WI./1e3,z_WI,TL_WI)
axis ij
shading flat
colormap(flipud(jet))
title('WI')
caxis([60 110])
xlabel('距离(km)')
ylabel('深度(m)')
c2 = colorbar;
ylabel(c2,'TL(dB)')

idx = 1:500;
mae_line = (mean(abs(TL_line_Ray(idx)-TL_line_WI(idx)).^1));

% idx_inf = 
mae_shd = (mean(abs(TL_Ray(:)-TL_WI(:)).^1));

%%
dataM = load('MunkB_CB_M.shd.mat','pressure','Pos');
TL_Ray_M = -20*log10(abs(squeeze(dataM.pressure)));
TL_Ray_M(1,:)=[];
TL_Ray_M(:,1)=[];

figure
pcolor(r_Ray./1e3,z_Ray,TL_Ray_M)
axis ij
shading flat
title('Ray-M')
colormap(flipud(jet))
caxis([60 110])
xlabel('距离(km)')
ylabel('深度(m)')
c3 = colorbar;
ylabel(c3,'TL(dB)')
%%
bellhopM_20260108_zjh('MunkB_CB_M')
%%
figure
plot(r_WI./1e3,TL_WI(1000,:),'k','linewidth',3)
axis ij
hold on
plot(r_Ray./1e3,TL_Ray(1000,:),'b','linewidth',2)
plot(r_Ray./1e3,TL_Ray_M(1000,:),'r','linewidth',2)
grid on
legend('WI','Bellhop','BellhopM')
xlabel('Range(km)')
ylabel('TL(dB)')