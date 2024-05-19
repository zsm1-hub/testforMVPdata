
function [I1,I2,J1,J2,time]=frame4320(lonpos1,lonpos2,latpos1,latpos2,TT,tindex,dot1,dot2,layer)

addpath(genpath('/sugon7/zsm/'));
addpath('/pacific3/SouthChinaSea2/');

DIRU=dir('/pacific3/SouthChinaSea2/U_py_960_1380_90_1000_T*.nc')
DIRV=dir('/pacific3/SouthChinaSea2/V_py_960_1380_90_1000_T*.nc')
DIRZ=dir('/pacific3/SouthChinaSea2/Eta_960_1380_10297.nc')
DIRT=dir('/pacific3/SouthChinaSea2/Theta_py_960_1380_90_1000_T*.nc')
DIRS=dir('/pacific3/SouthChinaSea2/Salt_py_960_1380_90_1000_T*.nc')

addpath('/pacific3/SouthChinaSea2/');
NX=960; NY=1380;
path0='/pacific3/SouthChinaSea2/grid/';
XC=readbin([path0,'XC_960x1380'],[NX NY]);
YC=readbin([path0,'YC_960x1380'],[NX NY]);
xx=XC';yy=YC';
% lonpos1=118;latpos1=23.5;
% lonpos2=120.5;latpos2=25;
% TT=4;tindex=120;
% layer=1;%1 is surface
% dot1=30;dot2=7;


dist=spheric_dist(latpos1,yy,lonpos1,xx);
[J11,I11]=find(dist==min(dist(:)));
I1=I11(1);J1=J11(1);

dist=spheric_dist(latpos2,yy,lonpos2,xx);
[J22,I22]=find(dist==min(dist(:)));
I2=I22(1);J2=J22(1);

fname1=DIRU(TT).name;
fname2=DIRV(TT).name;
fname3=DIRZ(1).name;
fname4=DIRT(TT).name;
fname5=DIRS(TT).name;

nc1=netcdf(fname1,'r');
nc2=netcdf(fname2,'r');
nc3=netcdf(fname3,'r');
nc4=netcdf(fname4,'r');
nc5=netcdf(fname5,'r');

if TT==1
    timecon=tindex;
else
    timecon=(TT-1)*1000+tindex;
end
time=datestr(datenum(2011,11,1)+timecon/24,'YYYY-mm-dd hh');

u=squeeze(nc1{'U'}(tindex,layer,:,:));
v=squeeze(nc2{'V'}(tindex,layer,:,:));
eta=squeeze(nc3{'Eta'}(timecon,:,:));
temp=squeeze(nc4{'Theta'}(tindex,layer,:,:));
salt=squeeze(nc5{'Salt'}(tindex,layer,:,:));

rho0=1025;T0=10;S0=35;TCOFF=1.7e-4;SCOEF=7.6e-4;
rho=rho0.*(1-TCOFF.*(temp-T0)+SCOEF.*(salt-S0));

u1=u(J1:J2,I1:I2);
v1=v(J1:J2,I1:I2);
salt1=salt(J1:J2,I1:I2);
temp1=temp(J1:J2,I1:I2);
eta1=eta(J1:J2,I1:I2);
rho1=rho(J1:J2,I1:I2);

xx1=xx(J1:J2,I1:I2);yy1=yy(J1:J2,I1:I2);

xframe=[xx(J1,I1) xx(J1,I2) xx(J2,I2) xx(J2,I1) xx(J1,I1)];
yframe=[yy(J1,I1) yy(J1,I2) yy(J2,I2) yy(J2,I1) yy(J1,I1)];

load mvppos.mat

target1=subplot(2,3,1)
pcolor(xx,yy,temp);shading flat;colorbar;hold on;
plot([xx(J1,I1) xx(J1,I2) xx(J2,I2) xx(J2,I1) xx(J1,I1)], [yy(J1,I1) yy(J1,I2) yy(J2,I2) yy(J2,I1) yy(J1,I1)],'k--','linewi',1.2);
line(lon0,lat0,'color','r','linewi',1.2)
line(lon,lat,'color','g','linewi',1,'linestyle',':')
quiver(xx(1:dot1:end,1:dot1:end),yy(1:dot1:end,1:dot1:end),...
    u(1:dot1:end,1:dot1:end),v(1:dot1:end,1:dot1:end),'w');
colortable = textread('temp1.txt');
%     colortable = flipud(colortable);
colormap(target1,colortable)
caxis([10 26]);
title('temp')



target2=subplot(2,3,4)
pcolor(xx1,yy1,temp1);shading flat;colorbar;hold on;
line(lon0,lat0,'color','r','linewi',1.2)
line(lon,lat,'color','g','linewi',1,'linestyle',':')
quiver(xx1(1:dot2:end,1:dot2:end),yy1(1:dot2:end,1:dot2:end),...
    u1(1:dot2:end,1:dot2:end),v1(1:dot2:end,1:dot2:end),'w');
colortable = textread('temp1.txt');
%     colortable = flipud(colortable);
colormap(target2,colortable)
caxis([10 26]);
title('temp')


target3=subplot(2,3,2)
pcolor(xx,yy,eta);shading flat;colorbar;hold on;
plot([xx(J1,I1) xx(J1,I2) xx(J2,I2) xx(J2,I1) xx(J1,I1)], [yy(J1,I1) yy(J1,I2) yy(J2,I2) yy(J2,I1) yy(J1,I1)],'k--','linewi',1.2);
line(lon0,lat0,'color','r','linewi',1.2)
line(lon,lat,'color','g','linewi',1,'linestyle',':')
quiver(xx(1:dot1:end,1:dot1:end),yy(1:dot1:end,1:dot1:end),...
    u(1:dot1:end,1:dot1:end),v(1:dot1:end,1:dot1:end),'w');
colortable = textread('CBR_wet.txt');
%     colortable = flipud(colortable);
colormap(target3,colortable)
caxis([-3 3]);
title('eta')


target4=subplot(2,3,5)
pcolor(xx1,yy1,eta1);shading flat;colorbar;hold on;
line(lon0,lat0,'color','r','linewi',1.2)
line(lon,lat,'color','g','linewi',1,'linestyle',':')
quiver(xx1(1:dot2:end,1:dot2:end),yy1(1:dot2:end,1:dot2:end),...
    u1(1:dot2:end,1:dot2:end),v1(1:dot2:end,1:dot2:end),'w');
colortable = textread('CBR_wet.txt');
%     colortable = flipud(colortable);
colormap(target4,colortable)
caxis([-3 3]);
title('eta')


target5=subplot(2,3,3)
pcolor(xx,yy,rho);shading flat;colorbar;hold on;
plot([xx(J1,I1) xx(J1,I2) xx(J2,I2) xx(J2,I1) xx(J1,I1)], [yy(J1,I1) yy(J1,I2) yy(J2,I2) yy(J2,I1) yy(J1,I1)],'k--','linewi',1.2);
line(lon0,lat0,'color','r','linewi',1.2)
line(lon,lat,'color','g','linewi',1,'linestyle',':')
quiver(xx(1:dot1:end,1:dot1:end),yy(1:dot1:end,1:dot1:end),...
    u(1:dot1:end,1:dot1:end),v(1:dot1:end,1:dot1:end),'w');
colortable = textread('matlab_jet.txt');
%     colortable = flipud(colortable);
colormap(target5,colortable)
caxis([1022 1024]);
title('rho')


target6=subplot(2,3,6)
pcolor(xx1,yy1,rho1);shading flat;colorbar;hold on;
line(lon0,lat0,'color','r','linewi',1.2)
line(lon,lat,'color','g','linewi',1,'linestyle',':')
quiver(xx1(1:dot2:end,1:dot2:end),yy1(1:dot2:end,1:dot2:end),...
    u1(1:dot2:end,1:dot2:end),v1(1:dot2:end,1:dot2:end),'w');
colortable = textread('matlab_jet.txt');
%     colortable = flipud(colortable);
colormap(target6,colortable)
caxis([1022 1024]);
title('rho')
suptitle([time])

return







