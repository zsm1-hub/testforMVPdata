clear all;close all;clc
fname='201903-台海数据汇交表(气象).xls'
addpath('D:\LIN2023\crocotools\COARE-algorithm-master\Matlab\COARE3.5');


B=readtable(fname);
(11+13/60)/24;
lat=table2array(B(2:end,6));
lon=table2array(B(2:end,7));

for ii=1:length(lat)
    lat1(ii)=str2double(lat{ii}(1:2))+str2double(lat{ii}(4:5))./60+...
        str2double(lat{ii}(7:13))./3600;
    lon1(ii)=str2double(lon{ii}(1:3))+str2double(lon{ii}(5:6))./60+...
        str2double(lon{ii}(8:14))./3600;
end

pressure=table2array(B(2:end,9));
temp=table2array(B(2:end,10));
windspeed=table2array(B(2:end,11));
winddir=table2array(B(2:end,12));
Relative_humidity=table2array(B(2:end,13));
hold on
load mvplonlat.mat
plot(lon_AB,lat_AB,'Color','b','LineWidth',2,'LineStyle','--')
plot(lon_CD,lat_CD,'Color','r','LineWidth',2,'LineStyle','--')
plot(lon1,lat1,'k');hold on


%% AB
% plot(lon1,lat1,'k');

dist1=spheric_dist(lat_AB(1),lat1,lon_AB(1),lon1);
I1=find(dist1==min(dist1))
dot=442;
lon_meteorology_AB=lon1(I1:I1+dot);
lat_meteorology_AB=lat1(I1:I1+dot);
plot(lon_meteorology_AB,lat_meteorology_AB,'k');
hold on;
plot(lon_AB,lat_AB,'Color','b','LineWidth',2,'LineStyle','--');

pressure_AB=pressure(I1:I1+dot);pressure_AB(pressure_AB>9000)=nan;
temp_AB=temp(I1:I1+dot);temp_AB(temp_AB>900)=nan;
windspeed_AB=windspeed(I1:I1+dot);windspeed_AB(windspeed_AB>90)=nan;
winddir_AB=winddir(I1:I1+dot);winddir_AB(winddir_AB>900)=nan;
Relative_humidity_AB=Relative_humidity(I1:I1+dot);
Relative_humidity_AB(Relative_humidity_AB>900)=nan;


windspd_AB=cell(1,size(lon_AB,2));
windDIR_AB=cell(1,size(lon_AB,2));
RH_AB=cell(1,size(lon_AB,2));
PRS_AB=cell(1,size(lon_AB,2));
TEMP_AB=cell(1,size(lon_AB,2));

for ii=1:length(windspd_AB)
    windspd_AB{ii}=NaN;
    windDIR_AB{ii}=NaN;
    RH_AB{ii}=NaN;
    PRS_AB{ii}=NaN;
    TEMP_AB{ii}=NaN;
end
for ii=1:size(lon_meteorology_AB,2)
    dist=spheric_dist(lat_meteorology_AB(ii),lat_AB,lon_meteorology_AB(ii),lon_AB);
    a=find(dist==min(dist));
    windspd_AB{a}=[windspd_AB{a} windspeed_AB(ii)];
    windDIR_AB{a}=[windDIR_AB{a} winddir_AB(ii)];
    RH_AB{a}=[RH_AB{a} Relative_humidity_AB(ii)];
    TEMP_AB{a}=[TEMP_AB{a} temp_AB(ii)];
    PRS_AB{a}=[PRS_AB{a} pressure_AB(ii)];
end

for ii=1:length(windspd_AB)
    windspd_AB1(ii)=nanmean(windspd_AB{ii});
    windDIR_AB1(ii)=nanmean(windDIR_AB{ii});
    RH_AB1(ii)=nanmean(RH_AB{ii});
    TEMP_AB1(ii)=nanmean(TEMP_AB{ii});
    PRS_AB1(ii)=nanmean(PRS_AB{ii});
end
windspd_AB2=fillmissing(windspd_AB1,'linear');
windDIR_AB2=fillmissing(windDIR_AB1,'linear');
RH_AB2=fillmissing(RH_AB1,'linear');
TEMP_AB2=fillmissing(TEMP_AB1,'linear');
PRS_AB2=fillmissing(PRS_AB1,'linear');
% 
% plot(lon_AB,windspd_AB2);
% hold on;
% plot(lon_AB,windspd_AB1,'r')

[uwind2,vwind2]=wind_vec(windspd_AB2,windDIR_AB2);
[uwind1,vwind1]=wind_vec(windspd_AB1,windDIR_AB1);

dot1=5;
hold on;
plot(lon_AB,lat_AB,'Color','k','LineWidth',2,'LineStyle','-');
quiver(lon_AB(1:dot1:end),lat_AB(1:dot1:end),uwind2(1:dot1:end),vwind2(1:dot1:end), ...
    1.2,'Color','k');

dot1=5;
hold on;
plot(lon_AB,lat_AB,'Color','r','LineWidth',2,'LineStyle','-');
quiver(lon_AB(1:dot1:end),lat_AB(1:dot1:end),uwind1(1:dot1:end),vwind1(1:dot1:end), ...
    1.2,'Color','k');


%% CD


dist1=spheric_dist(lat_CD(1),lat1,lon_CD(1),lon1);
I1=min(find(dist1==min(dist1)))

dist2=spheric_dist(lat_CD(end),lat1,lon_CD(end),lon1);
I2=min(find(dist2==min(dist2)))
lon_meteorology_CD=lon1(I1:I2);lat_meteorology_CD=lat1(I1:I2);
plot(lon_meteorology_CD,lat_meteorology_CD,'k');
hold on;
plot(lon_CD,lat_CD,'Color','r','LineWidth',2,'LineStyle','--');


pressure_CD=pressure(I1:I2);pressure_CD(pressure_CD>9000)=nan;
temp_CD=temp(I1:I2);temp_CD(temp_CD>900)=nan;
windspeed_CD=windspeed(I1:I2);windspeed_CD(windspeed_CD>90)=nan;
winddir_CD=winddir(I1:I2);winddir_CD(winddir_CD>900)=nan;
Relative_humidity_CD=Relative_humidity(I1:I2);
Relative_humidity_CD(Relative_humidity_CD>900)=nan;

windspd_CD=cell(1,size(lon_CD,2));
windDIR_CD=cell(1,size(lon_CD,2));
RH_CD=cell(1,size(lon_CD,2));
PRS_CD=cell(1,size(lon_CD,2));
TEMP_CD=cell(1,size(lon_CD,2));

for ii=1:length(windspd_CD)
    windspd_CD{ii}=NaN;
    windDIR_CD{ii}=NaN;
    RH_CD{ii}=NaN;
    PRS_CD{ii}=NaN;
    TEMP_CD{ii}=NaN;
end
for ii=1:size(lon_meteorology_CD,2)
    dist=spheric_dist(lat_meteorology_CD(ii),lat_CD,lon_meteorology_CD(ii),lon_CD);
    a=find(dist==min(dist));
    windspd_CD{a}=[windspd_CD{a} windspeed_CD(ii)];
    windDIR_CD{a}=[windDIR_CD{a} winddir_CD(ii)];
    RH_CD{a}=[RH_CD{a} Relative_humidity_CD(ii)];
    TEMP_CD{a}=[TEMP_CD{a} temp_CD(ii)];
    PRS_CD{a}=[PRS_CD{a} pressure_CD(ii)];
end

for ii=1:length(windspd_CD)
    windspd_CD1(ii)=nanmean(windspd_CD{ii});
    windDIR_CD1(ii)=nanmean(windDIR_CD{ii});
    RH_CD1(ii)=nanmean(RH_CD{ii});
    TEMP_CD1(ii)=nanmean(TEMP_CD{ii});
    PRS_CD1(ii)=nanmean(PRS_CD{ii});
end
windspd_CD2=fillmissing(windspd_CD1,'linear');
windDIR_CD2=fillmissing(windDIR_CD1,'linear');
RH_CD2=fillmissing(RH_CD1,'linear');
TEMP_CD2=fillmissing(TEMP_CD1,'linear');
PRS_CD2=fillmissing(PRS_CD1,'linear');

% 

plot(lon_CD,windspd_CD2);
hold on;
plot(lon_CD,windspd_CD1,'r');

[uwind2,vwind2]=wind_vec(windspd_CD2,windDIR_CD2);
[uwind1,vwind1]=wind_vec(windspd_CD1,windDIR_CD1);

dot1=5;
hold on;
plot(lon_CD,lat_CD,'Color','r','LineWidth',2,'LineStyle','-');
quiver(lon_CD(1:dot1:end),lat_CD(1:dot1:end),uwind2(1:dot1:end),vwind2(1:dot1:end), ...
    1.2,'Color','k');

dot1=5;
hold on;
plot(lon_CD,lat_CD,'Color','r','LineWidth',2,'LineStyle','-');
quiver(lon_CD(1:dot1:end),lat_CD(1:dot1:end),uwind1(1:dot1:end),vwind1(1:dot1:end), ...
    1.2,'Color','k');


[uwindalong,uwindcross,angle]=angle_mvp(lon_CD',lat_CD',uwind2',vwind2');
plot(lon_CD,uwindcross)

Cd=2.2e-3;


% plot(lon_CD,PRS_CD2)

% u=ws.*sin((wd+180)/360*2*pi);
% v=ws.*cos((wd+180)/360*2*pi);
% 
% % tmp=270.0-atan2(v,u)*180.0/pi;
% % dir=mod(tmp,360.0);
% % ws=sqrt(u.^2+v.^2);

