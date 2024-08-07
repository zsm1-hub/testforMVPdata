%%
clear all;close all;clc;
addpath(genpath('E:\ROMS学习\download_data_process\submeso\analysis\GSW\seawater\seawater'));
addpath('F:\TWS_Acrobat\TWS_Acrobat\TWS_Acrobat\')
addpath('E:\ROMS学习\download_data_process\submeso\initial')
addpath('E:\ROMS学习\download_data_process\submeso\analysis\taiwan')
% addpath('E:\ROMS学习\download_data_process\colorbar\colorbar_NCL\colorbar_NCL')
addpath('D:\colorbar\colorbar_NCL');
addpath('D:\LIN2023\getdata')
addpath('D:\mmap\m_map')
% load zsmmvp1.mat
% load zsmmvp2.mat
load ADCPzsm.mat
% load W116E124S22N30_20190310T2100Z.mat 
load W116E124S22N30_20190311T0000Z.mat
load mvplonlat.mat
lon_AB=fliplr(lon_AB);lat_AB=fliplr(lat_AB);
lon_CD=fliplr(lon_CD);lat_CD=fliplr(lat_CD);


J_range=[26:51];
I_range=[26:51];
[lon_r,lat_r]=meshgrid(lon,lat);
layer=1;
temp1=squeeze(temp(:,:,layer));salt1=squeeze(salt(:,:,layer));
us=squeeze(u(:,:,layer));vs=squeeze(v(:,:,layer));
dot=3;
lonr=lon_r';latr=lat_r';

%% large field
%% temp
% f1=subplot(2,2,1)
[lona,lonb,lata,latb]=m_range_zsm(lon_r,lat_r);
m_proj('mercator','lon',[lona lonb],'lat',[lata latb]);
m_pcolor(lonr,latr,temp1);shading interp;
c=colorbar;
set(get(c,'title'),'string','\circC')
% colortable=textread('thelix.txt');
colortable=textread('MPL_gnuplot.txt');
colormap(colortable);
m_line(lon_AB,lat_AB,'color','m','linewi',1.5);
m_line(lon_CD,lat_CD,'color','m','linewi',1.5);
m_gshhs_h('patch',[.7 .7 .7]);hold on;
m_quiver(lonr(1:dot:end,1:dot:end),latr(1:dot:end,1:dot:end),...
    us(1:dot:end,1:dot:end),vs(1:dot:end,1:dot:end),1.5,'color',[.7 .7 .7]);
m_line([lonr(J_range(1),I_range(1)),lonr(J_range(end),I_range(1)),...
    lonr(J_range(end),I_range(end)),lonr(J_range(1),I_range(end)),...
    lonr(J_range(1),I_range(1))],...
    [latr(J_range(1),I_range(1)),latr(J_range(end),I_range(1)),...
    latr(J_range(end),I_range(end)),latr(J_range(1),I_range(end)),...
    latr(J_range(1),I_range(1))],'linewi',1.5,'color','k')
m_grid('box','off','linestyle','none','fontsize',12,'FontWeight','b');
title('SST');
set(gca,'fontsize',12,'FontWeight','b')
saveas(gcf,'fig1_large_sst','png')

%% salt
% f2=subplot(2,2,2)
m_proj('mercator','lon',[lona lonb],'lat',[lata latb]);
m_pcolor(lonr,latr,salt1);shading interp;
c=colorbar;
set(get(c,'title'),'string','PSU')
colortable=textread('MPL_YlGnBu.txt');
colormap(colortable);
caxis([32 35]);
% load mvplonlat.mat
m_line(lon_AB,lat_AB,'color','m','linewi',1.5);
m_line(lon_CD,lat_CD,'color','m','linewi',1.5);
m_gshhs_h('patch',[.7 .7 .7]);hold on;
m_quiver(lonr(1:dot:end,1:dot:end),latr(1:dot:end,1:dot:end),...
    us(1:dot:end,1:dot:end),vs(1:dot:end,1:dot:end),1.5,'color',[.7 .7 .7]);
m_line([lonr(J_range(1),I_range(1)),lonr(J_range(end),I_range(1)),...
    lonr(J_range(end),I_range(end)),lonr(J_range(1),I_range(end)),...
    lonr(J_range(1),I_range(1))],...
    [latr(J_range(1),I_range(1)),latr(J_range(end),I_range(1)),...
    latr(J_range(end),I_range(end)),latr(J_range(1),I_range(end)),...
    latr(J_range(1),I_range(1))],'linewi',1.5,'color','k')
m_grid('box','off','linestyle','none','fontsize',12,'FontWeight','b');
title('SSS');
set(gca,'fontsize',12,'FontWeight','b')


% saveas(gcf,'fig1_large_sst','png')
% saveas(gcf,'fig1_large_sss','png')



% 
% clear
% load data_D_E_horizontal_interval_100m.mat
% lon_AB=DAT.lon_0;lat_AB=DAT.lat_0;
% load data_F_G_horizontal_interval_100m.mat
% lon_CD=DAT.lon_0;lat_CD=DAT.lat_0;
% save('mvplonlat.mat','lon_AB','lat_AB','lon_CD','lat_CD');

%% small field
dot2=1;
ssh2=ssh(J_range,I_range);
us2=us(J_range,I_range);vs2=vs(J_range,I_range);
lonr2=lonr(J_range,I_range);latr2=latr(J_range,I_range);
[lona,lonb,lata,latb]=m_range_zsm(lonr2,latr2);
m_proj('mercator','lon',[lona lonb],'lat',[lata latb]);
m_pcolor(lonr2,latr2,ssh2);shading interp;
c=colorbar;
set(get(c,'title'),'string','m')
colortable=textread('MPL_RdBu.txt');
colormap(flipud(colortable));
caxis([0.3 0.65]);
% load mvplonlat.mat
m_line(lon_AB,lat_AB,'color','m','linewi',1.5);
m_line(lon_CD,lat_CD,'color','m','linewi',1.5);
m_gshhs_h('patch',[.7 .7 .7]);hold on;
m_quiver(lonr2(1:dot2:end,1:dot2:end),latr2(1:dot2:end,1:dot2:end),...
    us2(1:dot2:end,1:dot2:end),vs2(1:dot2:end,1:dot2:end),1.5,'color','k');

% 
color_rho=flipud(flipud(textread('matlab_jet.txt')))
load zsmmvp1.mat;
for ii=1:size(rho,2)
    pos=min(find((~isnan(rho(:,ii)))==1));
    pos=11;
    if isempty(pos)==1
        rho_mvp_CD(:,ii)=NaN;
        temp_mvp_CD(:,ii)=NaN;
    else
        rho_mvp_CD(:,ii)=rho(pos,ii);
        temp_mvp_CD(:,ii)=temp(pos,ii);
    end
end

load zsmmvp2.mat;
for ii=1:size(rho,2)
    pos=min(find((~isnan(rho(:,ii)))==1));
    pos=11;
    if isempty(pos)==1
        rho_mvp_AB(:,ii)=NaN;
    else
        rho_mvp_AB(:,ii)=rho(pos,ii);
    end
end

rhomax=max([rho_mvp_AB rho_mvp_CD]);rhomin=min([rho_mvp_AB rho_mvp_CD]);
% rhomax=24.5;rhomin=23;
% rhomax=24.5;rhomin=15;
colorta=linspace(rhomin,rhomax,length(color_rho));

for ii=1:size(rho_mvp_CD,2)
      aa=abs(rho_mvp_CD(:,ii)-colorta);
      % aa=abs(temp_mvp_CD(:,ii)-colorta);
      if isnan(min(aa))==1
          poscolor(ii)=1;
      else
          poscolor(ii)=find(aa==min(aa));
      end
      m_scatter(lon_CD(ii),lat_CD(ii),7,color_rho(poscolor(ii),:),'filled');hold on
      % 
      disp(ii)
end

for ii=1:size(rho_mvp_AB,2)
      aa=abs(rho_mvp_AB(:,ii)-colorta);
      % aa=abs(temp_mvp_CD(:,ii)-colorta);
      if isnan(min(aa))==1
          poscolor(ii)=1;
      else
          poscolor(ii)=find(aa==min(aa));
      end
      m_scatter(lon_AB(ii),lat_AB(ii),7,color_rho(poscolor(ii),:),'filled');hold on
      % 
      disp(ii)
end

m_text(lon_AB(1,1)-0.6,lat_AB(1,1),'transect AB','fontsize',12,'FontWeight','b');
m_text(lon_CD(1,1)-0.6,lat_CD(1,1),'transect CD','fontsize',12,'FontWeight','b');


m_grid('box','off','linestyle','none','fontsize',12,'FontWeight','b','tickstyle','dm');
% title('SSS');
set(gca,'fontsize',12,'FontWeight','b')

saveas(gcf,'fig1_small_ssh','png')

%%%%%%%test strain

dudx=v2rho_2d((us2(2:end,:)-us2(1:end-1,:))./...
    spheric_dist(latr2(2:end,:),latr2(1:end-1,:),lonr2(2:end,:),lonr2(1:end-1,:)));
dvdx=v2rho_2d((vs2(2:end,:)-vs2(1:end-1,:))./...
    spheric_dist(latr2(2:end,:),latr2(1:end-1,:),lonr2(2:end,:),lonr2(1:end-1,:)));
dvdy=u2rho_2d((vs2(:,2:end)-vs2(:,1:end-1))./...
    spheric_dist(latr2(:,2:end),latr2(:,1:end-1),lonr2(:,2:end),lonr2(:,1:end-1)));
dudy=u2rho_2d((us2(:,2:end)-us2(:,1:end-1))./...
    spheric_dist(latr2(:,2:end),latr2(:,1:end-1),lonr2(:,2:end),lonr2(:,1:end-1)));
Snormal=dudx-dvdy;
Sshear=dvdx+dudy;
Strain=sqrt(Snormal.^2+Sshear.^2);
relative_vor=dvdx-dudy;
Ro=relative_vor./6e-5;

[lona,lonb,lata,latb]=m_range_zsm(lonr2,latr2);
m_proj('mercator','lon',[lona lonb],'lat',[lata latb]);
m_pcolor(lonr2,latr2,Strain);shading interp;hold on
m_quiver(lonr2(1:dot2:end,1:dot2:end),latr2(1:dot2:end,1:dot2:end),...
    us2(1:dot2:end,1:dot2:end),vs2(1:dot2:end,1:dot2:end),1.5,'color','k');
% pcolor(lonr2,latr2,Ro);shading interp
m_line(lon_AB,lat_AB,'color','m','linewi',1.5);
m_line(lon_CD,lat_CD,'color','m','linewi',1.5);
colortable=textread('MPL_YlGnBu.txt');
colormap(flipud(colortable));
c=colorbar
m_gshhs_h('patch',[.7 .7 .7]);
title('strain');
m_grid('box','off','linestyle','none','fontsize',12,'FontWeight','b','tickstyle','dm');
% title('SSS');


set(gca,'fontsize',12,'FontWeight','b')
saveas(gcf,'fig1_small_strain','png')

%% aviso
addpath(genpath('D:'))
fname='AVISO.nc';
ncdisp(fname,'/','full')
nc=netcdf(fname,'r');
tindex=2;
ugo=nc{'ugos'}(tindex,:,:);vgo=nc{'vgos'}(tindex,:,:);
sla=nc{'sla'}(tindex,:,:);
lona=nc{'longitude'}(:);
lata=nc{'latitude'}(:);
[lona1,lata1]=meshgrid(lona,lata);

load ADCPzsm.mat
load mvplonlat.mat
lon_AB=fliplr(lon_AB);lat_AB=fliplr(lat_AB);
lon_CD=fliplr(lon_CD);lat_CD=fliplr(lat_CD);
dot2=1;
[lona,lonb,lata,latb]=m_range_zsm(lona1,lata1);
m_proj('mercator','lon',[lona lonb],'lat',[lata latb]);
m_pcolor(lona1,lata1,sla);shading interp;
c=colorbar;
set(get(c,'title'),'string','m')
colortable=textread('MPL_RdBu.txt');
colormap(flipud(colortable));
caxis([-0.15 0.15]);
% load mvplonlat.mat
m_line(lon_AB,lat_AB,'color','m','linewi',1.5);
m_line(lon_CD,lat_CD,'color','m','linewi',1.5);
m_gshhs_h('patch',[.7 .7 .7]);hold on;
m_quiver(lona1(1:dot2:end,1:dot2:end),lata1(1:dot2:end,1:dot2:end),...
    ugo(1:dot2:end,1:dot2:end),vgo(1:dot2:end,1:dot2:end),1.5,'color','k');

% 
color_rho=flipud(flipud(textread('matlab_jet.txt')))
load zsmmvp1.mat;
for ii=1:size(rho,2)
    pos=min(find((~isnan(rho(:,ii)))==1));
    pos=11;
    if isempty(pos)==1
        rho_mvp_CD(:,ii)=NaN;
        temp_mvp_CD(:,ii)=NaN;
    else
        rho_mvp_CD(:,ii)=rho(pos,ii);
        temp_mvp_CD(:,ii)=temp(pos,ii);
    end
end

load zsmmvp2.mat;
for ii=1:size(rho,2)
    pos=min(find((~isnan(rho(:,ii)))==1));
    pos=11;
    if isempty(pos)==1
        rho_mvp_AB(:,ii)=NaN;
    else
        rho_mvp_AB(:,ii)=rho(pos,ii);
    end
end

rhomax=max([rho_mvp_AB rho_mvp_CD]);rhomin=min([rho_mvp_AB rho_mvp_CD]);
% rhomax=24.5;rhomin=23;
% rhomax=24.5;rhomin=15;
colorta=linspace(rhomin,rhomax,length(color_rho));

for ii=1:size(rho_mvp_CD,2)
      aa=abs(rho_mvp_CD(:,ii)-colorta);
      % aa=abs(temp_mvp_CD(:,ii)-colorta);
      if isnan(min(aa))==1
          poscolor(ii)=1;
      else
          poscolor(ii)=find(aa==min(aa));
      end
      m_scatter(lon_CD(ii),lat_CD(ii),7,color_rho(poscolor(ii),:),'filled');hold on
      % 
      disp(ii)
end

for ii=1:size(rho_mvp_AB,2)
      aa=abs(rho_mvp_AB(:,ii)-colorta);
      % aa=abs(temp_mvp_CD(:,ii)-colorta);
      if isnan(min(aa))==1
          poscolor(ii)=1;
      else
          poscolor(ii)=find(aa==min(aa));
      end
      m_scatter(lon_AB(ii),lat_AB(ii),7,color_rho(poscolor(ii),:),'filled');hold on
      % 
      disp(ii)
end

% m_text(lon_AB(1,1)-0.6,lat_AB(1,1),'transect AB','fontsize',12,'FontWeight','b');
% m_text(lon_CD(1,1)-0.6,lat_CD(1,1),'transect CD','fontsize',12,'FontWeight','b');


m_grid('box','off','linestyle','none','fontsize',12,'FontWeight','b','tickstyle','dm');
% title('SSS');
set(gca,'fontsize',12,'FontWeight','b')

saveas(gcf,'Fig1_aviso','png')
%% TS
clear all;close all;clc
load zsmmvp1.mat
layer=1:8;wid=1:size(temp,2);
Tmax=max(temp(:));Tmin=min(temp(:));
Smax=max(salt(:));Smin=min(salt(:));
rhop=rho+1000;
T=linspace(Tmin,Tmax,1000);S=linspace(Smin,Smax,1000);
[T1,S1]=meshgrid(T,S);
TT=temp(layer,wid);SS=salt(layer,wid);

% scatter(SS(:),TT(:),[],'filled','markerfaceColor',[.7 .7 .7])
alpha=1.7e-4;
beta=7.6e-4;
Tm=10;Sm=35;
rho_contour=1025.*[1-alpha.*(T1-Tm)+beta.*(S1-Sm)];

% scatter(salt(:),temp(:),[],'filled','markerfaceColor',[.7 .7 .7])
scatter(SS(:),TT(:),[],'filled','markerfaceColor',[.7 .7 .7])
% scatter(SS(:),TT(:),[],'filled','markerfaceColor','m')

hold on
load zsmmvp2.mat
scatter(salt(:),temp(:),[],'filled','markerfaceColor','m')

[C,h] =contour(S1,T1,rho_contour,'color','k','linestyle','-','linewi',1.2,'showtext','on')
clabel(C,h,'LabelSpacing',400);
xlabel('salinity PSU');
ylabel('temperature \circC');
title('T/S diagram')
set(gca,'fontsize',12,'FontWeight','b');


saveas(gcf,'fig1_TS_diagram','png')


%% TS 改

clear all;close all;clc
addpath(genpath('E:\ROMS学习\download_data_process\submeso\analysis\GSW\seawater\seawater'));
addpath('F:\TWS_Acrobat\TWS_Acrobat\TWS_Acrobat\')
addpath('E:\ROMS学习\download_data_process\submeso\initial')
addpath('E:\ROMS学习\download_data_process\submeso\analysis\taiwan')
% addpath('E:\ROMS学习\download_data_process\colorbar\colorbar_NCL\colorbar_NCL')
addpath('D:\colorbar\colorbar_NCL');
load zsmmvp1.mat
% load zsmmvp2.mat
load ADCPzsm.mat
load categories.mat

%%%%1是CD/FG，2是AB/DE

dx=abs(x(1,1)-x(1,2)).*1e3;
dz=abs(z(2,1)-z(1,1));
g=9.81;


xres=500;zres=-2;
xdot=abs(xres)./100;
zdot=abs(zres)./0.5;
clear temp1;clear salt1;clear rho1;
for ii=1:floor((size(temp,2)-1)/xdot)
    temp1(:,ii)=nanmean(temp(:,(ii-1)*xdot+1:(ii-1)*xdot+xdot),2);
    salt1(:,ii)=nanmean(salt(:,(ii-1)*xdot+1:(ii-1)*xdot+xdot),2);
    rho1(:,ii)=nanmean(rho(:,(ii-1)*xdot+1:(ii-1)*xdot+xdot),2);
end
clear temp2;clear salt2;clear rho2;
for ii=1:floor((size(temp1,1)-1)/zdot)
    temp2(ii,:)=nanmean(temp1((ii-1)*zdot+1:(ii-1)*zdot+zdot,:),1);
    salt2(ii,:)=nanmean(salt1((ii-1)*zdot+1:(ii-1)*zdot+zdot,:),1);
    rho2(ii,:)=nanmean(rho1((ii-1)*zdot+1:(ii-1)*zdot+zdot,:),1);
end

Tmax=max(temp2(:));Tmin=min(temp2(:));
Smax=max(salt2(:));Smin=min(salt2(:));
rhop=rho2+1000;
T=linspace(Tmin,Tmax,1000);S=linspace(Smin,Smax,1000);
[T1,S1]=meshgrid(T,S);
TT_S=temp2.*mask_s_front_CD;SS_S=salt2.*mask_s_front_CD;
TT_T=temp2.*mask_t_front_CD;SS_T=salt2.*mask_t_front_CD;
TT_M=temp2.*mask_transion_CD;SS_M=salt2.*mask_transion_CD;



%%%%% plot
alpha=1.7e-4;
beta=7.6e-4;
Tm=10;Sm=35;
rho_contour=1025.*[1-alpha.*(T1-Tm)+beta.*(S1-Sm)];


scatter(SS_S(:),TT_S(:),[],'filled','markerfaceColor','b','Marker','^',...
    'MarkerFaceAlpha', 0.3, 'MarkerEdgeAlpha', 0.2);
hold on
scatter(SS_T(:),TT_T(:),[],'filled','markerfaceColor','r','Marker','^',...
    'MarkerFaceAlpha', 0.3, 'MarkerEdgeAlpha', 0.2);
scatter(SS_M(:),TT_M(:),[],'filled','markerfaceColor','g','Marker','^',...
    'MarkerFaceAlpha', 0.3, 'MarkerEdgeAlpha', 0.2);


% hold on
% load zsmmvp2.mat
% scatter(salt(:),temp(:),[],'filled','markerfaceColor','m')

[C,h] =contour(S1,T1,rho_contour,'color','k','linestyle','-','linewi',1.2,'showtext','on')
clabel(C,h,'LabelSpacing',400);
xlabel('salinity PSU');
ylabel('temperature \circC');
title('T/S diagram')
set(gca,'fontsize',12,'FontWeight','b');


saveas(gcf,'fig1_TS_diagram5002','png')
