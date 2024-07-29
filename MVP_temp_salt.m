%%%%%%%%%%%%%%%
clear all;close all;clc
addpath(genpath('E:\ROMSÑ§Ï°\download_data_process\submeso\analysis\GSW\seawater\seawater'));
addpath('F:\TWS_Acrobat\TWS_Acrobat\TWS_Acrobat\')
addpath('E:\ROMSÑ§Ï°\download_data_process\submeso\initial')
addpath('E:\ROMSÑ§Ï°\download_data_process\colorbar\colorbar_NCL\colorbar_NCL')
% load('data_D_E_horizontal_interval_100m.mat');
% load('data_F_G_horizontal_interval_100m.mat');

load ADCP_AB_CD.mat
%%% AB is DE is zsmmvp2
%%% CD is FG is zsmmvp1
load('data_F_G_horizontal_interval_100m.mat');

%%

[Ds,De]=meshgrid(DAT.distance_0,-DAT.P_0);
Ds1=(Ds(1:end-1,1:end-1)+Ds(1:end-1,2:end))./2;
De1=(De(1:end-1,1:end-1)+De(2:end,1:end-1))./2;

PotDen=sw_dens(DAT.SA_0,DAT.pt_0,0);

ds_low=0:.1:max(DAT.distance_0);
de_low=-(0:.5:max(DAT.P_0));
[Ds_low,De_low]=meshgrid(ds_low,de_low);
Ds_low1=(Ds_low(1:end-1,1:end-1)+Ds_low(1:end-1,2:end))./2;
De_low1=(De_low(1:end-1,1:end-1)+De_low(2:end,1:end-1))./2;

PotDen_low=griddata(Ds,De,PotDen,Ds_low,De_low);
PotDen_low=fliplr(PotDen_low-1000);

PotTemp_low=griddata(Ds,De,DAT.pt_0,Ds_low,De_low);
PotTemp_low=fliplr(PotTemp_low);

PotSalt_low=griddata(Ds,De,DAT.SA_0,Ds_low,De_low);
PotSalt_low=fliplr(PotSalt_low);


Buo_low=(-9.8/1025).*(PotDen_low-1000);
N2_low=(Buo_low(2:end,:)-Buo_low(1:end-1,:))./...
    (De_low(2:end,:)-De_low(1:end-1,:));
N2_low=(N2_low(:,1:end-1)+N2_low(:,2:end))./2;
M2_low=(Buo_low(:,2:end)-Buo_low(:,1:end-1))./...
    (Ds_low(:,2:end)-Ds_low(:,1:end-1))./1000;
M2_low=(M2_low(1:end-1,:)+M2_low(2:end,:))./2;
f=sw_f(24.5);
Rib_low=((f.^2).*N2_low)./((M2_low).^2);
Rib_index=Rib_low;
Rib_index(Rib_index>=1)=1;
Rib_index(Rib_index>=0 & Rib_index<1)=2;
Rib_index(Rib_index>=-1 & Rib_index<0)=3;
Rib_index(Rib_index<-1)=4;



figure(1)
clf

f1=subplot(5,2,1)
pcolor(Ds_low,De_low,PotDen_low);
shading interp
hold on
contour(Ds_low,De_low,PotDen_low,[23.7,23.9,24],'k','LineWidth',2);
colorbar;
caxis([23.4,24.2]);
colormap("jet");
xticks(0:5:40);
yticks(-50:10:0);
ylabel('Depth (m)');
% A1=get(f1,'pos');
% A1(1)=A1(1)-0.05;
% A1(2)=A1(2)+0.035;
% A1(3)=A1(3)-0.00;
% A1(4)=A1(4)+0.03;
% set(f1,'pos',A1);
set(gca,'fontsize',12,'fontweight','bold');
title('density')

f2=subplot(5,2,3)
% pcolor(Ds_low,De_low,PotTemp_low);
contourf(Ds_low,De_low,PotTemp_low);
shading interp
hold on
contour(Ds_low,De_low,PotDen_low,[23.7,23.9,24],'k','LineWidth',2);
colorbar;
% caxis([23.4,24.2]);
colormap("jet");
xticks(0:5:40);
yticks(-50:10:0);
ylabel('Depth (m)');
% A1=get(f1,'pos');
% A1(1)=A1(1)-0.05;
% A1(2)=A1(2)+0.035;
% A1(3)=A1(3)-0.00;
% A1(4)=A1(4)+0.03;
% set(f1,'pos',A1);
set(gca,'fontsize',12,'fontweight','bold');
title('temp')

f3=subplot(5,2,5)
pcolor(Ds_low,De_low,PotSalt_low);
% contourf(Ds_low,De_low,PotSalt_low);
shading interp
hold on
contour(Ds_low,De_low,PotDen_low,[23.7,23.9,24],'k','LineWidth',2);
colorbar;
% caxis([23.4,24.2]);
colormap("jet");
xticks(0:5:40);
yticks(-50:10:0);
ylabel('Depth (m)');
caxis([34 34.8])
set(gca,'fontsize',12,'fontweight','bold');
title('salt')

subplot(5,2,6)
plot(Ds_low(1,:),PotSalt_low(8,:))

alpha=1.7e-4;
beta=7.6e-4;;
[R,Tu]=get_Turner(alpha,beta,PotTemp_low,PotSalt_low,Ds_low)

xx=0.5.*(Ds_low(:,2:end)+Ds_low(:,1:end-1));
zz=0.5.*(De_low(:,2:end)+De_low(:,1:end-1));
custom_colormap = [1 0 0; 1 1 0; 0 0 1; 0 1 0];  % ºì¡¢»Æ¡¢À¶¡¢ÂÌ

f4=subplot(5,2,7)
pcolor(xx,zz,Tu);
shading interp
hold on
contour(Ds_low,De_low,PotDen_low,[23.7,23.9,24],'k','LineWidth',2);
c=colorbar;
Tu(Tu==0)=nan;

caxis([-pi./2,pi./2]);
set(c,'ticks',[-pi./2,-pi/4,pi/4,pi/2],'ticklabels',{'-\pi/2','-\pi/4','\pi/4','\pi/2'});
colormap(f4,custom_colormap);

set(gca,'fontsize',12,'fontweight','bold');
title('Tu')

f5=subplot(5,2,9)
pcolor(xx,zz,R);
shading interp
hold on
contour(Ds_low,De_low,PotDen_low,[23.7,23.9,24],'k','LineWidth',2);
c=colorbar;
Tu(Tu==0)=nan;
colorbar;
colortable = textread('NCV_blue_red.txt');
colormap(f5,colortable)
% caxis([-pi./2,pi./2]);
% set(c,'ticks',[-pi./2,-pi/4,pi/4,pi/2],'ticklabels',{'-\pi/2','-\pi/4','\pi/4','\pi/2'});
set(gca,'fontsize',12,'fontweight','bold');
title('R')
caxis([-5 5])

subplot(5,2,10)
plot(xx(1,:),squeeze(Tu(30,:)),'b');hold on;
plot(xx(1,:),ones(1,size(R,2)).*pi/4,'k');
plot(xx(1,:),-ones(1,size(R,2)).*pi/4,'r');
title('vertical avg Tu')


f6=subplot(5,2,8)
N2=-9.8./1025.*v2rho_2d((PotDen_low(1:end-1,:)-PotDen_low(2:end,:))./(De_low(1:end-1,:)-De_low(2:end,:)));
N2_real=N2;
N2(N2<0)=1e-24;
% pcolor(Ds_low,De_low,log(N2));
pcolor(Ds_low,De_low,N2_real);
shading interp
hold on
contour(Ds_low,De_low,PotDen_low,[23.7,23.9,24],'k','LineWidth',2);
c=colorbar;
textname='temp_19lev.txt'
% 'GMT_ocean.txt'
colortable = textread(textname);
colormap(f6,colortable)
% caxis([-pi./2,pi./2]);
% set(c,'ticks',[-pi./2,-pi/4,pi/4,pi/2],'ticklabels',{'-\pi/2','-\pi/4','\pi/4','\pi/2'});
set(gca,'fontsize',12,'fontweight','bold');
title('N2')
% caxis([-20 -2])
caxis([-1e-4 1e-4])


f7=subplot(5,2,6)
bx=-9.8./1025.*u2rho_2d((PotDen_low(:,2:end)+1000-PotDen_low(:,1:end-1)-1000)./(Ds_low(:,2:end)-Ds_low(:,1:end-1)))./1e3;
pcolor(Ds_low,De_low,log(abs(bx)));
pcolor(Ds_low,De_low,(bx));
shading interp
hold on
contour(Ds_low,De_low,PotDen_low,[23.7,23.9,24],'k','LineWidth',2);
c=colorbar;
textname='matlab_jet.txt'
% 'GMT_ocean.txt'
% temp_19lev.txt
colortable = textread(textname);
colormap(f7,colortable)
% caxis([-pi./2,pi./2]);
% set(c,'ticks',[-pi./2,-pi/4,pi/4,pi/2],'ticklabels',{'-\pi/2','-\pi/4','\pi/4','\pi/2'});
set(gca,'fontsize',12,'fontweight','bold');
title('bx')
caxis([-2e-6 2e-6])
% caxis([-20 -10])

subplot(5,2,4)
Sx=-9.8./1025.*u2rho_2d((PotSalt_low(:,12:end)-PotSalt_low(:,1:end-11))./(Ds_low(:,12:end)-Ds_low(:,1:end-11)))./1e3;
Tx=-9.8./1025.*u2rho_2d((PotTemp_low(:,12:end)-PotTemp_low(:,1:end-11))./(Ds_low(:,12:end)-Ds_low(:,1:end-11)))./1e3;

plot(Sx(17,:),'b');hold on
plot(Tx(17,:),'r')

nanmean(Sx(17,208:end),2)
nanmean(Tx(17,1:140),2)

subplot(5,2,2)
plot(bx(17,:))


max(PotSalt_low(:))
min(PotSalt_low(:))
max(PotTemp_low(:))
min(PotTemp_low(:))

% saveas(gcf,'Tu','png');
z=De_low;x=Ds_low;
rho=PotDen_low;
temp=PotTemp_low;
salt=PotSalt_low;

% save('zsmmvp1.mat','alpha','beta','z','x','rho','temp','salt','f')
% save('zsmmvp2.mat','alpha','beta','z','x','rho','temp','salt','f')

%%% AB is DE is zsmmvp2
%%% CD is FG is zsmmvp1

lon=DAT.lon;lat=DAT.lat;
lon0=DAT.lon_0;lat0=DAT.lat_0;


[a,b]=meshgrid(lon,lat);
[lona,lonb,lata,latb]=m_range_zsm(a,b);


m_proj('mercator','lon',[117,122],'lat',[22,25]);
hold on;
m_line(lon,lat,'color','r','linewi',1.2);
m_line(lon0,lat0,'color','b')
m_gshhs_f('patch',[.7 .7 .7]);
% m_coast('patch',[.7 .7 .7]);
m_grid

% lon1=lon(2:end)-lon(1:end-1);
% lat1=lat(2:end)-lat(1:end-1);
% xxx=[0:250:40e3];
%for FG
xxx=[0:100:38.6e3];

u=DAT.u;v=DAT.v;
[ualong,uacross,angle]=angle_mvp(lon,lat,u,v);


[xddd2,u2]=avg_alongtrack(lon,lat,xxx,ualong);
[xddd2,v2]=avg_alongtrack(lon,lat,xxx,uacross);
lon1=lon;lat1=lat;
[xddd2,lon2]=avg_alongtrack(lon,lat,xxx,lon1);
[xddd2,lat2]=avg_alongtrack(lon,lat,xxx,lat1);

% plot(xxx(2:end),u2);hold on;
% plot(xxx(2:end),v2);

% u2=detrend(u2);
% v2=detrend(v2);
ualong1=fliplr(griddata(lon,lat,ualong,lon0,lat0));
uacross1=fliplr(griddata(lon,lat,uacross,lon0,lat0));


save('mvpposfg.mat','lon','lat','lon0','lat0','Ds_low','De_low','PotTemp_low','PotSalt_low','PotDen_low','u2','v2','xxx','lon2','lat2',...
    'u','v');

save('mvpadcpfg.mat','ualong1','uacross1');

% save('mvpposde.mat','lon','lat','lon0','lat0','Ds_low','De_low','PotTemp_low','PotSalt_low','PotDen_low','u2','v2','xxx','lon2','lat2',...
%     'u','v');
% 
% save('mvpadcpde.mat','ualong1','uacross1');

% PotDen=sw_pden(DAT.SA_0,DAT.pt_0,De,0);

PotDen=fliplr(PotDen);
Buo=(-9.8/1025).*PotDen;
N2=(Buo(2:end,:)-Buo(1:end-1,:))./...
    (De(2:end,:)-De(1:end-1,:));
N2=(N2(:,1:end-1)+N2(:,2:end))./2;
M2=(Buo(:,2:end)-Buo(:,1:end-1))./...
    (Ds(:,2:end)-Ds(:,1:end-1))./1000;
M2=(M2(1:end-1,:)+M2(2:end,:))./2;
f=sw_f(24.5);
Rib=((f.^2).*N2)./((M2).^2);

figure(2)
clf

subplot(2,2,1)
pcolor(Ds,-De,PotDen);
shading interp
hold on
contour(Ds,-De,PotDen,[1023.6,1023.9,1024],'k');
colorbar;
caxis([1023,1024.3]);
colormap("jet");

subplot(2,2,2)
% Rib(Rib<0)=1e-2;
% contourf(Ds1,-De1,N2,50,'linestyle','none');
% contourf(Ds1,-De1,log10(Rib),50,'linestyle','none');
pcolor(Ds1,-De1,N2);
shading interp
hold on
contour(Ds,-De,PotDen,[1023.6,1023.9,1024],'k');
colorbar;
caxis([-50e-3,50e-3]);
colormap("jet");

subplot(2,2,3)
% Rib(Rib<0)=1e-2;
% contourf(Ds1,-De1,N2,50,'linestyle','none');
% contourf(Ds1,-De1,log10(Rib),50,'linestyle','none');
pcolor(Ds1,-De1,M2);
shading interp
hold on
contour(Ds,-De,PotDen,[1023.6,1023.9,1024],'k');
colorbar;
caxis([-10e-2,10e-2]);
colormap("jet");

subplot(2,2,4)
Rib(Rib<0)=1e-2;
% contourf(Ds1,-De1,N2,50,'linestyle','none');
% contourf(Ds1,-De1,log10(Rib),50,'linestyle','none');
pcolor(Ds1,-De1,log10(1./Rib));
shading interp
hold on
contour(Ds,-De,PotDen,[1023.6,1023.9,1024],'k');
colorbar;
caxis([-3,3]);
colormap("jet");
%%

ds_low=0:.2:max(DAT.distance_0);
de_low=-(0:.5:max(DAT.P_0));
[Ds_low,De_low]=meshgrid(ds_low,de_low);
Ds_low1=(Ds_low(1:end-1,1:end-1)+Ds_low(1:end-1,2:end))./2;
De_low1=(De_low(1:end-1,1:end-1)+De_low(2:end,1:end-1))./2;
PotDen_low=griddata(Ds,De,PotDen,Ds_low,De_low);
PotDen_low=fliplr(PotDen_low-1000);
Buo_low=(-9.8/1025).*(PotDen_low-1000);
N2_low=(Buo_low(2:end,:)-Buo_low(1:end-1,:))./...
    (De_low(2:end,:)-De_low(1:end-1,:));
N2_low=(N2_low(:,1:end-1)+N2_low(:,2:end))./2;
M2_low=(Buo_low(:,2:end)-Buo_low(:,1:end-1))./...
    (Ds_low(:,2:end)-Ds_low(:,1:end-1))./1000;
M2_low=(M2_low(1:end-1,:)+M2_low(2:end,:))./2;
f=sw_f(24.5);
Rib_low=((f.^2).*N2_low)./((M2_low).^2);
Rib_index=Rib_low;
Rib_index(Rib_index>=1)=1;
Rib_index(Rib_index>=0 & Rib_index<1)=2;
Rib_index(Rib_index>=-1 & Rib_index<0)=3;
Rib_index(Rib_index<-1)=4;

figure(2)
clf

f1=subplot(2,2,1)
pcolor(Ds_low,De_low,PotDen_low);
shading interp
hold on
contour(Ds_low,De_low,PotDen_low,[23.6,23.9,24],'k');
colorbar;
caxis([23.5,24.2]);
colormap("jet");

f2=subplot(2,2,2)
pcolor(Ds_low1,De_low1,N2_low);
shading interp
hold on
contour(Ds_low,De_low,PotDen_low,[23.6,23.9,24],'k');
colorbar;
caxis([-1e-3,1e-3]);
colormap(f2,GMT_polar);

f3=subplot(2,2,3)
pcolor(Ds_low1,De_low1,M2_low);
shading interp
hold on
contour(Ds_low,De_low,PotDen_low,[23.6,23.9,24],'k');
colorbar;
caxis([-2e-6,2e-6]);
colormap(f3,GMT_polar);

f3=subplot(2,2,3)
pcolor(Ds_low1,De_low1,Rib_low);
shading interp
hold on
contour(Ds_low,De_low,PotDen_low,[23.7,23.9,24],'b','LineWidth',2);
% contour(Ds_low1,De_low1,Rib_low,[.25,1],'g','LineWidth',2);
colorbar;
caxis([-2,4]);
% colormap(f4,flipud(WhiteBlue));
colormap(f3,hsv);
% colormap(f3,(cool));
% colormap(f4,dia_color);

f4=subplot(2,2,4)

map=[1,0,0;1,1,0;0,1,0;0,0,1];
scatter(reshape(Ds_low1,[],1),reshape(De_low1,[],1),10,reshape(Rib_index,[],1),'s','filled');
c=colorbar;
caxis([1,4]);
colormap(f4,flipud(map));
set(c,'ticks',1:4,'ticklabels',{'BCI&stable','SI','GI&SI','GI'});
hold on
contour(Ds_low,De_low,PotDen_low,[23.7,23.9,24],'k','LineWidth',2);
box on

%%

figure(1)
clf

subplot(2,1,1)
pcolor(Ds,-De,fliplr(DAT.pt_0));
shading interp
hold on
contour(Ds,-De,fliplr(DAT.pt_0),[17,20],'k');
caxis([15,22]);
colormap("jet");

subplot(2,1,2)
pcolor(Ds,-De,fliplr(DAT.SA_0));
shading interp
hold on
contour(Ds,-De,fliplr(DAT.SA_0),[33,34],'k');
caxis([31,35]);
colormap("jet");

figure(2)
clf

scatter(DAT.lon,DAT.lat);
