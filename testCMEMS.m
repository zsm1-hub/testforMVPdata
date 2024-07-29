clear all;close all;clc
% fname='cmems_mod_glo_phy_my_0.083deg_P1D-m_1716049535443.nc';
fname='CMEMS2.nc';

ncdisp(fname,'/','full')
[lon,lat,uu,vv,zos]=read_CMEMS(fname)
[lon1,lat1]=meshgrid(lon,lat);
layer=1;dot1=3;

figure(1)
subplot(2,2,1)
contourf(lon1,lat1,zos);colorbar;hold on;
% caxis([])
quiver(lon1(1:dot1:end,1:dot1:end),lat1(1:dot1:end,1:dot1:end)...
    ,squeeze(uu(layer,1:dot1:end,1:dot1:end)),squeeze(vv(layer,1:dot1:end,1:dot1:end)),'k')
load mvppos.mat
line(lon0,lat0,'color','r','linewi',1.2)
title('Reanalysis','fontweight','b')
set(gca,'fontweight','b')

[u0c]=interp2(lon1,lat1,squeeze(uu(layer,:,:)),lon2,lat2);
[v0c]=interp2(lon1,lat1,squeeze(vv(layer,:,:)),lon2,lat2);

subplot(2,2,3)
quiver(lon2,lat2,u2,v2,'k');hold on;
quiver(lon2,lat2,u0c,v0c,'b');
R = corr2([u2;v2], [u0c;v0c])
% text(119,24.3,['R=',num2str(R)],'fontweight','b')
text(119.4,24.65,['R=',num2str(R)],'fontweight','b')
set(gca,'fontweight','b')





avisofname='AVISO2.nc'
ncdisp(avisofname,'/','full')

[lona,lata,ugos,vgos,sla]=read_AVISO(avisofname);
[lona1,lata1]=meshgrid(lona,lata);
dot2=1;

subplot(2,2,2)
contourf(lona1,lata1,sla);colorbar;hold on;
% caxis([])
quiver(lona1(1:dot2:end,1:dot2:end),lata1(1:dot2:end,1:dot2:end)...
    ,squeeze(ugos(1:dot2:end,1:dot2:end)),squeeze(vgos(1:dot2:end,1:dot2:end)),'k')
load mvppos.mat
line(lon0,lat0,'color','r','linewi',1.2)
title('AVISO','fontweight','b')
set(gca,'fontweight','b')


[u0a]=interp2(lona1,lata1,ugos,lon2,lat2);
[v0a]=interp2(lona1,lata1,vgos,lon2,lat2);

subplot(2,2,4)
quiver(lon2,lat2,u2,v2,'k');hold on;
quiver(lon2,lat2,u0a,v0a,'b');
R = corr2([u2;v2], [u0a;v0a])
% text(119,24.3,['R=',num2str(R)],'fontweight','b')
text(119.4,24.65,['R=',num2str(R)],'fontweight','b')

set(gca,'fontweight','b')
saveas(gcf,'Related_AVISO_CMEMS2','png')