clc;clear all;
freqList = [10, 12.5, 16, 20, 25, 31.5, 40, 50, 63, 80, ...
    100, 125, 160, 200, 250, 315, 400, 500, 630, 800, ...
    1000, 1250, 1600, 2000, 2500, 3150, 4000, 5000, 6300, 8000, 10000];
for num = 1:length(freqList)
    FreqDir = ['/workspace/zjh/Result/' num2str(freqList(num)) 'Hz/'];
    if exist(FreqDir, 'dir') == 0
        mkdir(FreqDir);
        disp(['创建文件夹: ' FreqDir]);
    else
        disp(['文件夹: ' FreqDir '存在'])
    end
end
tic
for ifreq = 15:length(freqList)
FreqDir = ['/workspace/zjh/Result/' num2str(freqList(ifreq)) 'Hz/'];
f0 = freqList(ifreq);   % 声源频率 Hz
sd = 1000;   % 声源深度 m
rhob = 1.8; % 海底密度 g/cm^3
cb = 1600;  % 海底声速 m/s
attnb = 0.8;% 海底衰减系数 单位需根据后面的ssp_option(3)设定，目前ssp_option(3)='W'，单位为dB/波长
rmax = 100; % 最远接收距离 km
rr = 0:0.1:rmax; % 距离方向接收水听器位置设置
rd = 0:1:5000;   % 深度方向接收水听器位置设置
% 声速剖面的深度向量zw
zw = [0 ;200;250;400;600;800;1000;1200;1400;1600;1800;2000;2200;2400;2600;2800;3000;3200;3400;3600;3800;4000;4200;4400;4600;4800;5000];
% 声速剖面的声速向量cw
cw = [1548.52 ;1530.29;1526.69;1517.78;1509.49;1504.30;1501.38;1500.14;1500.12;1501.02;1502.57;1504.62;1507.02;1509.69;1512.55;1515.56;1518.67;1521.85;1525.10;1528.38;1531.7;1535.04;1538.39;1541.76;1545.14;1548.52;1551.91];
zmax = 5000; % 最大接收深度
%% Write Env——bellhop程序输入环境文件设置和写入(后缀为.env文件)
envname = 'MunkB_CB_M'; % env文件名，可自定义
nmedia = 1; % 介质分层数，在bellhop中，这个参数一直是1
ssp_option = 'CVW'; %  第1个字母是SSP的插值类型选择，第2个字母是海面状态类型选择，第3个字母是衰减系数单位选择，具体详见《海洋声传播模拟（一）射线模拟.pdf》第9页
seabed_option = 'A'; % 海底状态类型选择，A表示流体半空间
sspbottom = [max(zw) cb 0 rhob attnb 0];  % 海底声学参数设置
runtype_option = 'CB'; % 声场计算类型选择。若算声场，第1字母为C表示相干声场，第1个字母若为I表示非相干声场；若算声线，只需要runtype_option = 'R';
nbeams = 0; % 声线数选择，设为0的话程序会自动计算需要的声线根数
alpha=[-60 60]; % 声线出射范围，这里设置的是-90~90
zbox = 5500;%ceil(zmax*1.1); % 声线追踪区域深度方向，通常需要比最深接收水听器深一些
rbox = 101;%ceil(rmax*1.1); % 声线追踪区域距离方向，通常需要比最远接收水听器远一些

fid = fopen([FreqDir envname  '.env'],'w');
fprintf(fid, '''%s.env'' \r\n',envname);
fprintf(fid, '%8.3f  !Freq \r\n',f0);
fprintf(fid, '%d \r\n',nmedia);
fprintf(fid, '''%s'' \r\n',ssp_option);
fprintf(fid, '0 0.0  %8.3f \r\n', max(zw));
m=length(zw(:));
for im=1:m
    fprintf(fid, '%8.3f %8.3f / \r\n',zw(im),cw(im));
end
fprintf(fid, '''%s'' 0.0 \r\n', seabed_option);
fprintf(fid, '%8.3f %8.3f %8.3f %8.3f %8.3f %8.3f \r\n',sspbottom(1),sspbottom(2),sspbottom(3),sspbottom(4),sspbottom(5),sspbottom(6));
fprintf(fid,'%g  \r\n',length(sd));
fprintf(fid,'%8.3f %8.3f / \r\n',sd(1),sd(end));
fprintf(fid,'%g  \r\n',length(rd));
fprintf(fid,'%8.3f %8.3f / \r\n',rd(1),rd(end));
fprintf(fid,'%g  \r\n',length(rr));
fprintf(fid,'%8.3f %8.3f / \r\n',rr(1),rr(end));
fprintf(fid,'''%s''                 !Run type:Ray/Coh/Inc/Sem \r\n',runtype_option);
fprintf(fid,'%d  0              !NBEAMS IBEAM \r\n',nbeams);
fprintf(fid,'%8.3f %8.3f /      !Alpha1,2 \r\n',alpha(1),alpha(end));
fprintf(fid,'0.0 %8.3f %8.3f \r\n',zbox,rbox);
fclose(fid);
% clear fid;
%% Run——运行
bellhopM_20260108_zjh([FreqDir envname]);
end
toc

