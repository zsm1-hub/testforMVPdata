clear all;close all;clc
addpath(genpath('E:\QMDownload\SoftMgr\LeapFTP3.0.1.46/'));
addpath('E:\ROMSѧϰ\download_data_process\colorbar\colorbar_NCL\colorbar_NCL');
fnamechl1='H08_20190310_2000_1H_ROC021_FLDK.02401_02401.nc'
% fnamechl2='H08_20190311_0000_1H_ROC021_FLDK.02401_02401.nc'
fnamesst='20190310000000-JAXA-L3C_GHRSST-SSTskin-H08_AHI-v1.2_daily-v02.0-fv01.0.nc';

ncdisp(fnamesst,'/','full')
nct=netcdf(fnamesst,'r');
sst=nct{'sea_surface_temperature'}(:);
spd=nct{'wind_speed'}(:);

add_offset=273.15;
scale_factor=0.01;
lon=nct{'lon'}(:);
lat=nct{'lat'}(:);
[lon1,lat1]=meshgrid(lon,lat);
lona=118
lonb=123
lata=22
latb=28
[J1,I1]=find(spheric_dist(lata,lat1,lona,lon1)==min(min(spheric_dist(lata,lat1,lona,lon1))));
[J2,I2]=find(spheric_dist(latb,lat1,lonb,lon1)==min(min(spheric_dist(latb,lat1,lonb,lon1))));

lon22=lon1(J2:J1,I1:I2);lat22=lat1(J2:J1,I1:I2);sst2=sst(J2:J1,I1:I2);spd2=spd(J2:J1,I1:I2);

sst3=sst2.*scale_factor;sst3(abs(sst3)>100)=nan;
target1=subplot(221)
[lonaa,lonbb,lataa,latbb]=m_range_zsm(lon22,lat22);
m_proj('mercator','lon',[lonaa lonbb],'lat',[lataa latbb]);
m_pcolor(lon22,lat22,sst3);shading interp;colorbar;
load mvppos.mat;hold on;
m_plot(lon2,lat2,'color','k','linewi',1.2)
colorbar;
colortable = textread('matlab_jet.txt');
colormap(target1,colortable);
m_gshhs_h('patch',[.7 .7 .7]);
m_grid;


target2=subplot(222)
[lonaa,lonbb,lataa,latbb]=m_range_zsm(lon22,lat22);
m_proj('mercator','lon',[lonaa lonbb],'lat',[lataa latbb]);
m_pcolor(lon22,lat22,spd2);shading interp;colorbar;
load mvppos.mat;hold on;
m_plot(lon2,lat2,'color','k','linewi',1.2)
colorbar;
colortable = textread('matlab_jet.txt');
colormap(target2,jet(11));
m_gshhs_h('patch',[.7 .7 .7]);
m_grid;


ncdisp(fnamechl1,'/','full')
ncchl1=netcdf(fnamechl1,'r');
lonc=ncchl1{'longitude'}(:);
latc=ncchl1{'latitude'}(:);
chlor_a=ncchl1{'chlor_a'}(:).*0.004;% mg/m^3
chlor_a(abs(chlor_a)>130)=nan;

[lonc1,latc1]=meshgrid(lonc,latc);
lona=118
lonb=123
lata=22
latb=28
[Jc1,Ic1]=find(spheric_dist(lata,latc1,lona,lonc1)==min(min(spheric_dist(lata,latc1,lona,lonc1))));
[Jc2,Ic2]=find(spheric_dist(latb,latc1,lonb,lonc1)==min(min(spheric_dist(latb,latc1,lonb,lonc1))));
loncc=lonc1(Jc2:Jc1,Ic1:Ic2);latcc=latc1(Jc2:Jc1,Ic1:Ic2);chlor_a2=chlor_a(Jc2:Jc1,Ic1:Ic2);

target3=subplot(223)
[lonaa,lonbb,lataa,latbb]=m_range_zsm(loncc,latcc);
m_proj('mercator','lon',[lonaa lonbb],'lat',[lataa latbb]);
m_pcolor(loncc,latcc,log(chlor_a2));shading interp;colorbar;
load mvppos.mat;hold on;
m_plot(lon2,lat2,'color','k','linewi',1.2)
colorbar;
colortable = textread('MPL_winter.txt');
colormap(target3,(colortable));
m_gshhs_h('patch',[.7 .7 .7]);
m_grid;
caxis([-2 4])


