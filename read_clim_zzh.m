clear all;close all;clc
addpath(genpath('/sugon7/zsm/'));
addpath('/leader/user/zzh/oldroms_nestingv2/res/clim/clim3Y/');
fname='roms_twschinav2clim3Y_his.n_0003.nc';
tindex=10;
[u,v,w,rho,temp,salt,zeta,h]=get_param_zsm(fname,tindex);
nc=netcdf(fname,'r');
lon_rho=nc{'lon_rho'}(:);
lat_rho=nc{'lat_rho'}(:);
chlorophyll=nc{'chlorophyll'}(tindex,:,:,:);
sustr=u2rho_2d(nc{'sustr'}(tindex,:,:));
svstr=v2rho_2d(nc{'svstr'}(tindex,:,:));

chlorophyll(isnan(temp))=nan;
sustr(isnan(squeeze(temp(end,:,:))))=nan;
svstr(isnan(squeeze(temp(end,:,:))))=nan;
sustr(abs(sustr)>100)=nan;svstr(abs(svstr)>100)=nan;


% [~,~,~,rho]=get_hslice_ParamCroco(fname,z,tindex);
rho=linear_EOS(temp,salt);



%% temp
textname='temp1.txt'
vmin=0;vmax=25;varname='temp'
subplot(321)
pcolor(lon_rho,lat_rho,squeeze(temp(end,:,:)));shading interp;hold on;
%         contour(x,y,rhoh,'linewi',1.2,'color','g','linestyle','-');
quiver(lon_rho(1:dot2:end,1:dot1:end),lat_rho(1:dot2:end,1:dot1:end),...
sustr(1:dot2:end,1:dot1:end),svstr(1:dot2:end,1:dot1:end),0.6);
colorbar;
colortable = textread(textname);
% eval(['colormap(','target',num2str(pos1),',colortable)'])
colormap(colortable)
% caxis([-2e-4 2e-4]);
caxis([vmin,vmax]);
title([varname])
load mvpposde.mat
line(lon0,lat0,'linewi',1.5,'color','k')
text(lon0(1),lat0(1),'E')
text(lon0(end),lat0(end),'D')
load mvpposfg.mat
line(lon0,lat0,'linewi',1.5,'color','k')
text(lon0(1),lat0(1),'G')
text(lon0(end),lat0(end),'F')


subplot(323)
load mvpposde.mat
pcolor(Ds_low,De_low,PotTemp_low);shading interp;colorbar
caxis([16 22])
title('model section D-E')


subplot(324)
load mvpposde.mat
lonsec=[lon0(end) lon0(1)];latsec=[lat0(end) lat0(1)];
[X,Z,temp1]=get_section(fname,fname,lonsec,latsec,'temp',tindex);
pcolor(X,Z,temp1);shading interp;colorbar
caxis([12.6 14.1])
title('mvp section D-E')

subplot(325)
load mvpposfg.mat
pcolor(Ds_low,De_low,PotTemp_low);shading interp;colorbar
caxis([12 23])
title('model section F-G')

subplot(326)
load mvpposfg.mat
lonsec=[lon0(end) lon0(1)];latsec=[lat0(end) lat0(1)];
[X,Z,temp1]=get_section(fname,fname,lonsec,latsec,'temp',tindex);
pcolor(X,Z,temp1);shading interp;colorbar
caxis([12.6 14.1])
title('section F-G')

saveas(gcf,'clim_mvp_temp','png')


%% rho
textname='matlab_jet.txt'
vmin=18;vmax=26;varname='rho'
subplot(321)
pcolor(lon_rho,lat_rho,squeeze(rho(end,:,:))-1000);shading interp;hold on;
%         contour(x,y,rhoh,'linewi',1.2,'color','g','linestyle','-');
quiver(lon_rho(1:dot2:end,1:dot1:end),lat_rho(1:dot2:end,1:dot1:end),...
sustr(1:dot2:end,1:dot1:end),svstr(1:dot2:end,1:dot1:end),0.6);
colorbar;
colortable = textread(textname);
% eval(['colormap(','target',num2str(pos1),',colortable)'])
colormap(colortable)
% caxis([-2e-4 2e-4]);
caxis([vmin,vmax]);
title([varname])
load mvpposde.mat
line(lon0,lat0,'linewi',1.5,'color','k')
text(lon0(1),lat0(1),'E')
text(lon0(end),lat0(end),'D')

load mvpposfg.mat
line(lon0,lat0,'linewi',1.5,'color','k')
text(lon0(1),lat0(1),'G')
text(lon0(end),lat0(end),'F')

subplot(323)
load mvpposde.mat
pcolor(Ds_low,De_low,PotDen_low);shading interp;colorbar
caxis([23.6 24.2])
title('model section D-E')


subplot(324)
load mvpposde.mat
lonsec=[lon0(end) lon0(1)];latsec=[lat0(end) lat0(1)];
[X,Z,rho1]=get_section(fname,fname,lonsec,latsec,'rho',tindex);
pcolor(X,Z,rho1);shading interp;colorbar
caxis([23.8 26.5])
title('mvp section D-E')

subplot(325)
load mvpposfg.mat
pcolor(Ds_low,De_low,PotDen_low);shading interp;colorbar
caxis([22.5 24.2])
title('model section F-G')

subplot(326)
load mvpposfg.mat
lonsec=[lon0(end) lon0(1)];latsec=[lat0(end) lat0(1)];
[X,Z,rho1]=get_section(fname,fname,lonsec,latsec,'rho',tindex);
pcolor(X,Z,rho1);shading interp;colorbar
caxis([24 26])
title('section F-G')
saveas(gcf,'clim_mvp_rho','png')
%% bio
% CHLOROPHYLL
textname='MPL_BuGn.txt';dot1=55;dot2=dot1*2;
vmin=0;vmax=3;varname='chlorophyll'
target1=subplot(321)
pcolor(lon_rho,lat_rho,log(squeeze(chlorophyll(end,:,:))));
shading interp;hold on;
contour(lon_rho,lat_rho,squeeze(rho(end,:,:))-1000,'color','k');
quiver(lon_rho(1:dot2:end,1:dot1:end),lat_rho(1:dot2:end,1:dot1:end),...
sustr(1:dot2:end,1:dot1:end),svstr(1:dot2:end,1:dot1:end),0.6);
%         contour(x,y,rhoh,'linewi',1.2,'color','g','linestyle','-');
colorbar;
colortable = textread(textname);
% eval(['colormap(','target',num2str(pos1),',colortable)'])
colormap(target1,colortable)
% caxis([-2e-4 2e-4]);
caxis([vmin,vmax]);
title([varname])
load mvpposde.mat
line(lon0,lat0,'linewi',1.5,'color','k')
text(lon0(1),lat0(1),'E')
text(lon0(end),lat0(end),'D')
load mvpposfg.mat
line(lon0,lat0,'linewi',1.5,'color','k')
text(lon0(1),lat0(1),'G')
text(lon0(end),lat0(end),'F')



target3=subplot(323)
load mvpposde.mat
lonsec=[lon0(end) lon0(1)];latsec=[lat0(end) lat0(1)];
[X,Z,chlorophyll1]=get_section(fname,fname,lonsec,latsec,'chlorophyll',tindex);
[X,Z,rho1]=get_section(fname,fname,lonsec,latsec,'rho',tindex);
pcolor(X,Z,chlorophyll1);shading interp;colorbar;hold on;
contour(X,Z,rho1,'linewi',1.5,'color','k');
colortable = textread('MPL_BuGn.txt');
colormap(target3,colortable)
caxis([0 3])
title('model section D-E')

target5=subplot(325)
load mvpposde.mat
lonsec=[lon0(end) lon0(1)];latsec=[lat0(end) lat0(1)];
[X,Z, no31]=get_section(fname,fname,lonsec,latsec,'NO3',tindex);
pcolor(X,Z,no31);shading interp;colorbar;hold on;
contour(X,Z,rho1,'linewi',1.5,'color','k');
colortable = textread('MPL_cool.txt');
colormap(target5,colortable)
caxis([0 1.5])
title('model section D-E')

target4=subplot(324)
load mvpposfg.mat
lonsec=[lon0(end) lon0(1)];latsec=[lat0(end) lat0(1)];
[X,Z,chlorophyll1]=get_section(fname,fname,lonsec,latsec,'chlorophyll',tindex);
[X,Z,rho1]=get_section(fname,fname,lonsec,latsec,'rho',tindex);
pcolor(X,Z,chlorophyll1);shading interp;colorbar;hold on;
contour(X,Z,rho1,'linewi',1.5,'color','k');
colortable = textread('MPL_BuGn.txt');
colormap(target4,colortable)
caxis([0 3])
title('model section F-G')

target6=subplot(326)
load mvpposfg.mat
lonsec=[lon0(end) lon0(1)];latsec=[lat0(end) lat0(1)];
[X,Z, no31]=get_section(fname,fname,lonsec,latsec,'NO3',tindex);
pcolor(X,Z,no31);shading interp;colorbar;hold on;
contour(X,Z,rho1,'linewi',1.5,'color','k');
colortable = textread('MPL_cool.txt');
colormap(target6,colortable)
caxis([0 1.5])
title('model section F-G')
saveas(gcf,'clim_mvp_bio','png')
