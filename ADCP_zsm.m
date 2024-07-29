%%%%%%%%%%%%%%%
clear all;close all;clc
addpath(genpath('E:\ROMS学习\download_data_process\submeso\analysis\GSW\seawater\seawater'));
addpath('F:\TWS_Acrobat\TWS_Acrobat\TWS_Acrobat\')
addpath('E:\ROMS学习\download_data_process\submeso\initial')
addpath('E:\ROMS学习\download_data_process\colorbar\colorbar_NCL\colorbar_NCL')
% load('data_D_E_horizontal_interval_100m.mat');
% load('data_F_G_horizontal_interval_100m.mat');

load ADCP_AB_CD.mat
%%% AB is DE is zsmmvp2
%%% CD is FG is zsmmvp1

%% CD---FG
load('data_F_G_horizontal_interval_100m.mat');

plot(lon_CD,lat_CD,'color','r');hold on;
plot(DAT.lon_0,DAT.lat_0,'Color','k')

%% 插值
for i=1:size(u_CD,2)
    vel1=squeeze(u_CD(:,i));vel2=squeeze(v_CD(:,i));
    u_CD1(:,i)=griddata(lon_CD,lat_CD,vel1,DAT.lon_0,DAT.lat_0,"nearest");
    v_CD1(:,i)=griddata(lon_CD,lat_CD,vel2,DAT.lon_0,DAT.lat_0,"nearest");
    disp(i)
end
[Ds,De]=meshgrid(DAT.distance_0,-DAT.P_0);
zz=-depth;xx=Ds(1,:);
[xx1,zz1]=meshgrid(xx,zz);
u_CD2=(griddata(xx1,zz1,u_CD1',Ds,De));
v_CD2=(griddata(xx1,zz1,v_CD1',Ds,De));

%% PLOT TEST
subplot(2,2,1)
pcolor(Ds,De,(u_CD2));shading interp;
colortable=textread('NCV_blue_red.txt');
colormap(colortable);colorbar;
caxis([-1 1]);
title('u_CD','interpreter','none')

subplot(2,2,2)
pcolor(Ds,De,(v_CD2));shading interp;
colortable=textread('NCV_blue_red.txt');
colormap(colortable);colorbar;
caxis([-1 1]);
title('v_CD','interpreter','none')


%% 速度投影
[ualong_CD,uacross_CD,angle]=Vel_project(DAT.lon_0',DAT.lat_0',u_CD2',v_CD2');
ualong_CD1=fliplr(ualong_CD');
uacross_CD1=fliplr(uacross_CD');
%% plot test
subplot(2,2,1)
pcolor(Ds,De,ualong_CD1);shading interp;
colortable=textread('NCV_blue_red.txt');
colormap(colortable);colorbar;
caxis([-1 1]);
title('ualong_CD','interpreter','none')

subplot(2,2,2)
pcolor(Ds,De,uacross_CD1);shading interp;
colortable=textread('NCV_blue_red.txt');
colormap(colortable);colorbar;
caxis([-1 1]);
title('uacross_CD','interpreter','none')
%% AB---CD
load('data_D_E_horizontal_interval_100m.mat');

plot(lon_AB,lat_AB,'color','r');hold on;
plot(DAT.lon_0,DAT.lat_0,'Color','k')

%% 插值
for i=1:size(u_AB,2)
    vel1=squeeze(u_AB(:,i));vel2=squeeze(v_AB(:,i));
    u_AB1(:,i)=griddata(lon_AB,lat_AB,vel1,DAT.lon_0,DAT.lat_0,"nearest");
    v_AB1(:,i)=griddata(lon_AB,lat_AB,vel2,DAT.lon_0,DAT.lat_0,"nearest");
    disp(i)
end
[Ds,De]=meshgrid(DAT.distance_0,-DAT.P_0);
zz=-depth;xx=Ds(1,:);
[xx1,zz1]=meshgrid(xx,zz);
u_AB2=griddata(xx1,zz1,u_AB1',Ds,De);
v_AB2=griddata(xx1,zz1,v_AB1',Ds,De);
%% 速度投影
[ualong_AB,uacross_AB,angle]=Vel_project(DAT.lon_0',DAT.lat_0',u_AB2',v_AB2');
ualong_AB1=fliplr(ualong_AB');
uacross_AB1=fliplr(uacross_AB');
%% plot test
subplot(2,2,3)
pcolor(Ds,De,ualong_AB1);shading interp;
colortable=textread('NCV_blue_red.txt');
colormap(colortable);colorbar;
caxis([-1 1]);
title('ualong_AB','interpreter','none')

subplot(2,2,4)
pcolor(Ds,De,uacross_AB1);shading interp;
colortable=textread('NCV_blue_red.txt');
colormap(colortable);colorbar;
caxis([-1 1]);
title('uacross_AB','interpreter','none')


save('ADCPzsm.mat','ualong_AB1','uacross_AB1','ualong_CD1','uacross_CD1')