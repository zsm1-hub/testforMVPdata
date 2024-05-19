clear all;clc;close all;
U=[];V=[];
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

scrsize=get(0,'screensize');%% 是为了获得屏幕大小，Screensize是一个4元素向量[left,bottom, width, height]
scrsize(3)=scrsize(3)./2; %
set(gcf,'position',scrsize);% 用获得的screensize向量设置figure的position属性，实现最大化的目的

TT=4;
% tindex=120;

for tindex=96:144
    lonpos1=118;latpos1=23.5;
    lonpos2=120.5;latpos2=25;
    layer=1;%1 is surface
    dot1=30;dot2=7;
%     [I1,I2,J1,J2,time]=frame4320(lonpos1,lonpos2,latpos1,latpos2,TT,tindex,dot1,dot2,layer)
    [Ds_low,De_low,PotDen_low,time]=vertical_interpolate(lonpos1,lonpos2,latpos1,latpos2,TT,tindex);
    saveas(gcf,['LLC4320_MVP_vertical',time],'png')
    clf
    %%%%%%%Ds为x,De为z


end
% suptitle(['Days = ',num2str(tindot*6/24,'%2.2f')]);
        
