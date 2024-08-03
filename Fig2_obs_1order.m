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



% dbdx=u2rho_2d((rho(:,2:end)-rho(:,1:end-1))./dx)./1025.*-g;
% dtdx=u2rho_2d((temp(:,2:end)-temp(:,1:end-1))./dx);
% dsdx=u2rho_2d((salt(:,2:end)-salt(:,1:end-1))./dx);
% 
% effect_rhox_t=-alpha.*dtdx;
% effect_rhox_s=beta.*dsdx;
% [mld]=get_mld_taiwanPIO(rho,z);



%% 粗化 CD
clf;
load zsmmvp1.mat
xres=500;zres=-2;
xdot=abs(xres)./100;
zdot=abs(zres)./0.5;
clear temp1;clear salt1;clear rho1;
clear uacross1; clear ualong1;
for ii=1:floor((size(temp,2)-1)/xdot)
    temp1(:,ii)=nanmean(temp(:,(ii-1)*xdot+1:(ii-1)*xdot+xdot),2);
    salt1(:,ii)=nanmean(salt(:,(ii-1)*xdot+1:(ii-1)*xdot+xdot),2);
    rho1(:,ii)=nanmean(rho(:,(ii-1)*xdot+1:(ii-1)*xdot+xdot),2);
        uacross1(:,ii)=nanmean(uacross_CD1(:,(ii-1)*xdot+1:(ii-1)*xdot+xdot),2);
        ualong1(:,ii)=nanmean(ualong_CD1(:,(ii-1)*xdot+1:(ii-1)*xdot+xdot),2);
    % uacross1(:,ii)=nanmean(uacross_AB1(:,(ii-1)*xdot+1:(ii-1)*xdot+xdot),2);
    % ualong1(:,ii)=nanmean(ualong_AB1(:,(ii-1)*xdot+1:(ii-1)*xdot+xdot),2);
end
clear temp2;clear salt2;clear rho2;
clear uacross2;clear ualong2;
for ii=1:floor((size(temp1,1)-1)/zdot)
    temp2(ii,:)=nanmean(temp1((ii-1)*zdot+1:(ii-1)*zdot+zdot,:),1);
    salt2(ii,:)=nanmean(salt1((ii-1)*zdot+1:(ii-1)*zdot+zdot,:),1);
    rho2(ii,:)=nanmean(rho1((ii-1)*zdot+1:(ii-1)*zdot+zdot,:),1);
    uacross2(ii,:)=nanmean(uacross1((ii-1)*zdot+1:(ii-1)*zdot+zdot,:),1);
    ualong2(ii,:)=nanmean(ualong1((ii-1)*zdot+1:(ii-1)*zdot+zdot,:),1);
end

x1=0:xres:size(temp2,2).*xres-xres;
z1=[0:zres:size(temp2,1).*zres-zres]';
[x2,z2]=meshgrid(x1,z1);
x2=x2./1e3;

g=9.8
dbdx2=u2rho_2d((rho2(:,2:end)-rho2(:,1:end-1))./xres)./1025.*-g;
dtdx2=u2rho_2d((temp2(:,2:end)-temp2(:,1:end-1))./xres);
dsdx2=u2rho_2d((salt2(:,2:end)-salt2(:,1:end-1))./xres);
dbdz2=v2rho_2d(((rho2(1:end-1,:)-rho2(2:end,:))./abs(zres))*-g./1025);
dtdz2=v2rho_2d((temp2(1:end-1,:)-temp2(2:end,:))./abs(zres));
dsdz2=v2rho_2d((salt2(1:end-1,:)-salt2(2:end,:))./abs(zres));
dbdz2(dbdz2<4e-10)=4e-10;


pycnal=0.1;

figure;
left=0.15;
bot=0.8;
width=0.8;
height=0.17;
zpos=0.18;
colorcon='k'
load mld.mat
mld=mld_CD;

mask1=abs(rho2-24.1)<pycnal/1;
%%补充
mask1(1:3,41:44)=1;


% 创建第一个坐标系
f1=axes('Position', [left, bot, width, height]); % 左下角坐标 (0.1, 0.7)，宽度和高度各为 0.25
pcolor(x2,z2,temp2);shading interp;colorbar;hold on;
contour(x2,z2,rho2,[22.5:pycnal:24],'linewi',.5,'linestyle','-','color',colorcon);
contour(x2,z2,mask1,[1 1],'linewi',1.5,'linestyle','--','color',colorcon);
plot(mld_CDx(1,:),mld,'color','r','linestyle','--','LineWidth',1.5);
colortable=textread('MPL_gnuplot.txt');
colormap(f1,colortable);
c=colorbar;
set(c,'ytick',[15 18 21])
caxis([13 24])
ylabel('depth [m]');
text(1.5,-50,'T')
text(17,5,'Transect CD','FontWeight','b')
text(5,-20,'S front','FontWeight','b')
text(12,-35,'transion zone','FontWeight','b')
text(30,-20,'T front','FontWeight','b')
% text(17,5,'Transect AB','FontWeight','b')

set(gca,'xtick',[],'ytick',[-50 -30 -10])
set(gca,'fontsize',10,'fontweight','b');

f2=axes('Position', [left, bot-zpos*1, width, height]); 
pcolor(x2,z2,salt2);shading interp;colorbar;hold on;
contour(x2,z2,rho2,[22.5:pycnal:24],'linewi',.5,'linestyle','-','color',colorcon);
contour(x2,z2,mask1,[1 1],'linewi',1.5,'linestyle','--','color',colorcon);
plot(mld_CDx(1,:),mld,'color','r','linestyle','--','LineWidth',1.5);
colortable=textread('MPL_YlGnBu.txt');
colormap(f2,(colortable));
caxis([30.4 34.8])
ylabel('depth [m]');
text(1.5,-50,'S')
c=colorbar;
% set(c,'ytick',[31 32.2 33.4])
set(gca,'xtick',[],'ytick',[-50 -30 -10])
set(gca,'fontsize',10,'fontweight','b');


f3=axes('Position', [left, bot-zpos*2, width, height]); 
pcolor(x2,z2,rho2);shading interp;colorbar;hold on;
contour(x2,z2,rho2,[22.5:pycnal:24],'linewi',.5,'linestyle','-','color',colorcon);
contour(x2,z2,mask1,[1 1],'linewi',1.5,'linestyle','--','color',colorcon);
plot(mld_CDx(1,:),mld,'color','r','linestyle','--','LineWidth',1.5);
colortable=textread('matlab_jet.txt');
colormap(f3,colortable);
% caxis([22.5 24])
caxis([23 24.5])
c=colorbar;
set(c,'ytick',[23.5 24])
ylabel('depth [m]');
text(1.5,-50,'\sigma')
set(gca,'xtick',[],'ytick',[-50 -30 -10])
set(gca,'fontsize',10,'fontweight','b');


f4=axes('Position', [left, bot-zpos*3, width, height]); 
pcolor(x2,z2,uacross2);shading interp;colorbar;hold on;
contour(x2,z2,rho2,[22.5:pycnal:24],'linewi',.5,'linestyle','-','color',colorcon);
contour(x2,z2,mask1,[1 1],'linewi',1.5,'linestyle','--','color',colorcon);
plot(mld_CDx(1,:),mld,'color','r','linestyle','--','LineWidth',1.5);
colortable=textread('MPL_RdBu.txt');
colormap(f4,flipud(colortable));
caxis([-.5 .5])
ylabel('depth [m]');
text(1.5,-50,'Vel_cross','interpreter','none')
set(gca,'xtick',[],'ytick',[-50 -30 -10])
set(gca,'fontsize',10,'fontweight','b');


f5=axes('Position', [left, bot-zpos*4, width, height]); 
pcolor(x2,z2,ualong2);shading interp;colorbar;hold on;
contour(x2,z2,rho2,[22.5:pycnal:24],'linewi',.5,'linestyle','-','color',colorcon);
contour(x2,z2,mask1,[1 1],'linewi',1.5,'linestyle','--','color',colorcon);
plot(mld_CDx(1,:),mld,'color','r','linestyle','--','LineWidth',1.5);
colortable=textread('MPL_RdBu.txt');
colormap(f5,flipud(colortable));
caxis([-.5 .5])
h = xlabel('km');
set(h, 'Position', [20, -65, -0.5], 'Units', 'normalized');
ylabel('depth [m]');
text(1.5,-50,'Vel_along','interpreter','none')
set(gca,'ytick',[-50 -30 -10],'xtick',[0 10 20 30 40])
set(gca,'fontsize',10,'fontweight','b');


saveas(gcf,'grid5002CD','png')
% saveas(gcf,'grid5001AB','png')
%% CD原始
%% 原始一阶
figure;
left=0.15;
bot=0.8;
width=0.8;
height=0.17;
zpos=0.18;

% 创建第一个坐标系
load mld.mat;
mld=mld_CD;
f1=axes('Position', [left, bot, width, height]); % 左下角坐标 (0.1, 0.7)，宽度和高度各为 0.25
pcolor(x,z,temp);shading interp;colorbar;hold on;
contour(x,z,rho,[22.5:0.1:24],'linewi',.5,'linestyle','-','color','k');
contour(x2,z2,mask1,[1 1],'linewi',1.5,'linestyle','--','color',colorcon);
plot(mld_CDx(1,:),mld,'color','r','linestyle','--','LineWidth',1.5);
colortable=textread('MPL_gnuplot.txt');
colormap(f1,colortable);
c=colorbar;
set(c,'ytick',[15 18 21])
caxis([13 23])
ylabel('depth [m]');
text(1.5,-50,'T')
text(17,5,'Transect CD','FontWeight','b')
text(5,-20,'S front','FontWeight','b')
text(12,-35,'transion zone','FontWeight','b')
text(30,-20,'T front','FontWeight','b')
set(gca,'xtick',[],'ytick',[-50 -30 -10])
set(gca,'fontsize',10,'fontweight','b');

f2=axes('Position', [left, bot-zpos*1, width, height]); 
pcolor(x,z,salt);shading interp;colorbar;hold on;
contour(x,z,rho,[22.5:0.1:24],'linewi',.5,'linestyle','-','color','k');
contour(x2,z2,mask1,[1 1],'linewi',1.5,'linestyle','--','color',colorcon);

plot(mld_CDx(1,:),mld,'color','r','linestyle','--','LineWidth',1.5);
colortable=textread('MPL_YlGnBu.txt');
colormap(f2,flipud(colortable));
caxis([30 34.8])
ylabel('depth [m]');
text(1.5,-50,'S')
c=colorbar;
% set(c,'ytick',[31 32.2 33.4])
set(gca,'xtick',[],'ytick',[-50 -30 -10])
set(gca,'fontsize',10,'fontweight','b');


f3=axes('Position', [left, bot-zpos*2, width, height]); 
pcolor(x,z,rho);shading interp;colorbar;hold on;
contour(x,z,rho,[22.5:0.1:24],'linewi',.5,'linestyle','-','color','k');
contour(x2,z2,mask1,[1 1],'linewi',1.5,'linestyle','--','color',colorcon);

plot(mld_CDx(1,:),mld,'color','r','linestyle','--','LineWidth',1.5);
colortable=textread('matlab_jet.txt');
colormap(f3,colortable);
% caxis([22.5 24])
caxis([23 24.5])
c=colorbar;
set(c,'ytick',[23.5 24])
ylabel('depth [m]');
text(1.5,-50,'\sigma')
set(gca,'xtick',[],'ytick',[-50 -30 -10])
set(gca,'fontsize',10,'fontweight','b');


f4=axes('Position', [left, bot-zpos*3, width, height]); 
pcolor(x,z,uacross_CD1);shading interp;colorbar;hold on;
contour(x,z,rho,[22.5:0.1:24],'linewi',.5,'linestyle','-','color','k');
contour(x2,z2,mask1,[1 1],'linewi',1.5,'linestyle','--','color',colorcon);
plot(mld_CDx(1,:),mld,'color','r','linestyle','--','LineWidth',1.5);
colortable=textread('MPL_RdBu.txt');
colormap(f4,flipud(colortable));
caxis([-.5 .5])
ylabel('depth [m]');
text(1.5,-50,'Vel_cross','interpreter','none')
set(gca,'xtick',[],'ytick',[-50 -30 -10])
set(gca,'fontsize',10,'fontweight','b');


f5=axes('Position', [left, bot-zpos*4, width, height]); 
pcolor(x,z,ualong_CD1);shading interp;colorbar;hold on;
contour(x,z,rho,[22.5:0.1:24],'linewi',.5,'linestyle','-','color','k');
contour(x2,z2,mask1,[1 1],'linewi',1.5,'linestyle','--','color',colorcon);
plot(mld_CDx(1,:),mld,'color','r','linestyle','--','LineWidth',1.5);
colortable=textread('MPL_RdBu.txt');
colormap(f5,flipud(colortable));
caxis([-.5 .5])
h = xlabel('km');
set(h, 'Position', [20, -65, -0.5], 'Units', 'normalized');
ylabel('depth [m]');
text(1.5,-50,'Vel_along','interpreter','none')
set(gca,'ytick',[-50 -30 -10],'xtick',[0 10 20 30 40])
set(gca,'fontsize',10,'fontweight','b');


saveas(gcf,'grid100CD','png')
% saveas(gcf,'grid100AB','png')


%% 粗化 AB
clf
load zsmmvp2.mat
xres=500;zres=-2;
xdot=abs(xres)./100;
zdot=abs(zres)./0.5;
clear temp1;clear salt1;clear rho1;
clear uacross1; clear ualong1;
for ii=1:floor((size(temp,2)-1)/xdot)
    temp1(:,ii)=nanmean(temp(:,(ii-1)*xdot+1:(ii-1)*xdot+xdot),2);
    salt1(:,ii)=nanmean(salt(:,(ii-1)*xdot+1:(ii-1)*xdot+xdot),2);
    rho1(:,ii)=nanmean(rho(:,(ii-1)*xdot+1:(ii-1)*xdot+xdot),2);
        % uacross1(:,ii)=nanmean(uacross_CD1(:,(ii-1)*xdot+1:(ii-1)*xdot+xdot),2);
        % ualong1(:,ii)=nanmean(ualong_CD1(:,(ii-1)*xdot+1:(ii-1)*xdot+xdot),2);
    uacross1(:,ii)=nanmean(uacross_AB1(:,(ii-1)*xdot+1:(ii-1)*xdot+xdot),2);
    ualong1(:,ii)=nanmean(ualong_AB1(:,(ii-1)*xdot+1:(ii-1)*xdot+xdot),2);
end

clear temp2;clear salt2;clear rho2;
clear uacross2;clear ualong2;
for ii=1:floor((size(temp1,1)-1)/zdot)
    temp2(ii,:)=nanmean(temp1((ii-1)*zdot+1:(ii-1)*zdot+zdot,:),1);
    salt2(ii,:)=nanmean(salt1((ii-1)*zdot+1:(ii-1)*zdot+zdot,:),1);
    rho2(ii,:)=nanmean(rho1((ii-1)*zdot+1:(ii-1)*zdot+zdot,:),1);
    uacross2(ii,:)=nanmean(uacross1((ii-1)*zdot+1:(ii-1)*zdot+zdot,:),1);
    ualong2(ii,:)=nanmean(ualong1((ii-1)*zdot+1:(ii-1)*zdot+zdot,:),1);
end

x1=0:xres:size(temp2,2).*xres-xres;
z1=[0:zres:size(temp2,1).*zres-zres]';
[x2,z2]=meshgrid(x1,z1);
x2=x2./1e3;

g=9.8
dbdx2=u2rho_2d((rho2(:,2:end)-rho2(:,1:end-1))./xres)./1025.*-g;
dtdx2=u2rho_2d((temp2(:,2:end)-temp2(:,1:end-1))./xres);
dsdx2=u2rho_2d((salt2(:,2:end)-salt2(:,1:end-1))./xres);
dbdz2=v2rho_2d(((rho2(1:end-1,:)-rho2(2:end,:))./abs(zres))*-g./1025);
dtdz2=v2rho_2d((temp2(1:end-1,:)-temp2(2:end,:))./abs(zres));
dsdz2=v2rho_2d((salt2(1:end-1,:)-salt2(2:end,:))./abs(zres));
dbdz2(dbdz2<4e-10)=4e-10;


pycnal=0.1;

figure;
left=0.15;
bot=0.8;
width=0.8;
height=0.17;
zpos=0.18;
colorcon='k'
load mld.mat
mld=mld_AB;

mask1=abs(rho2-24.2)<pycnal/1;
%%补充
mask1(1:2,52:63)=1;
% 创建第一个坐标系
f1=axes('Position', [left, bot, width, height]); % 左下角坐标 (0.1, 0.7)，宽度和高度各为 0.25
pcolor(x2,z2,temp2);shading interp;colorbar;hold on;
contour(x2,z2,rho2,[22.5:pycnal:24],'linewi',.5,'linestyle','-','color',colorcon);
contour(x2,z2,mask1,[1 1],'linewi',1.5,'linestyle','--','color',colorcon);
plot(mld_ABx(1,:),mld,'color','r','linestyle','--','LineWidth',1.5);
colortable=textread('MPL_gnuplot.txt');
colormap(f1,colortable);
c=colorbar;
set(c,'ytick',[15 18 21])
caxis([13 24])
ylabel('depth [m]');
text(1.5,-50,'T')
text(5,-20,'S front','FontWeight','b')
text(21,-35,'transion zone','FontWeight','b')
text(35,-20,'T front','FontWeight','b')
% text(17,5,'Transect CD','FontWeight','b')
text(17,5,'Transect AB','FontWeight','b')

set(gca,'xtick',[],'ytick',[-50 -30 -10])
set(gca,'fontsize',10,'fontweight','b');

f2=axes('Position', [left, bot-zpos*1, width, height]); 
pcolor(x2,z2,salt2);shading interp;colorbar;hold on;
contour(x2,z2,rho2,[22.5:pycnal:24],'linewi',.5,'linestyle','-','color',colorcon);
contour(x2,z2,mask1,[1 1],'linewi',1.5,'linestyle','--','color',colorcon);
plot(mld_ABx(1,:),mld,'color','r','linestyle','--','LineWidth',1.5);
colortable=textread('MPL_YlGnBu.txt');
colormap(f2,(colortable));
caxis([30.4 34.8])
ylabel('depth [m]');
text(1.5,-50,'S')
c=colorbar;
% set(c,'ytick',[31 32.2 33.4])
set(gca,'xtick',[],'ytick',[-50 -30 -10])
set(gca,'fontsize',10,'fontweight','b');


f3=axes('Position', [left, bot-zpos*2, width, height]); 
pcolor(x2,z2,rho2);shading interp;colorbar;hold on;
contour(x2,z2,rho2,[22.5:pycnal:24],'linewi',.5,'linestyle','-','color',colorcon);
contour(x2,z2,mask1,[1 1],'linewi',1.5,'linestyle','--','color',colorcon);
plot(mld_ABx(1,:),mld,'color','r','linestyle','--','LineWidth',1.5);
colortable=textread('matlab_jet.txt');
colormap(f3,colortable);
% caxis([22.5 24])
caxis([23 24.5])
c=colorbar;
set(c,'ytick',[23.5 24])
ylabel('depth [m]');
text(1.5,-50,'\sigma')
set(gca,'xtick',[],'ytick',[-50 -30 -10])
set(gca,'fontsize',10,'fontweight','b');


f4=axes('Position', [left, bot-zpos*3, width, height]); 
pcolor(x2,z2,uacross2);shading interp;colorbar;hold on;
contour(x2,z2,rho2,[22.5:pycnal:24],'linewi',.5,'linestyle','-','color',colorcon);
contour(x2,z2,mask1,[1 1],'linewi',1.5,'linestyle','--','color',colorcon);
plot(mld_ABx(1,:),mld,'color','r','linestyle','--','LineWidth',1.5);
colortable=textread('MPL_RdBu.txt');
colormap(f4,flipud(colortable));
caxis([-.5 .5])
ylabel('depth [m]');
text(1.5,-50,'Vel_cross','interpreter','none')
set(gca,'xtick',[],'ytick',[-50 -30 -10])
set(gca,'fontsize',10,'fontweight','b');


f5=axes('Position', [left, bot-zpos*4, width, height]); 
pcolor(x2,z2,ualong2);shading interp;colorbar;hold on;
contour(x2,z2,rho2,[22.5:pycnal:24],'linewi',.5,'linestyle','-','color',colorcon);
contour(x2,z2,mask1,[1 1],'linewi',1.5,'linestyle','--','color',colorcon);
plot(mld_ABx(1,:),mld,'color','r','linestyle','--','LineWidth',1.5);
colortable=textread('MPL_RdBu.txt');
colormap(f5,flipud(colortable));
caxis([-.5 .5])
h = xlabel('km');
set(h, 'Position', [20, -65, -0.5], 'Units', 'normalized');
ylabel('depth [m]');
text(1.5,-50,'Vel_along','interpreter','none')
set(gca,'ytick',[-50 -30 -10],'xtick',[0 10 20 30 40])
set(gca,'fontsize',10,'fontweight','b');


saveas(gcf,'grid5002AB','png')

%% AB原始
clf
load zsmmvp2.mat
figure;
left=0.15;
bot=0.8;
width=0.8;
height=0.17;
zpos=0.18;
mld=mld_AB;
mask1=abs(rho2-24.2)<pycnal/1;
%%补充
mask1(1:2,52:63)=1;
% 创建第一个坐标系
f1=axes('Position', [left, bot, width, height]); % 左下角坐标 (0.1, 0.7)，宽度和高度各为 0.25
pcolor(x,z,temp);shading interp;colorbar;hold on;
contour(x,z,rho,[22.5:0.15:24],'linewi',.5,'linestyle','-','color','k');
contour(x2,z2,mask1,[1 1],'linewi',1.5,'linestyle','--','color',colorcon);
plot(mld_ABx(1,:),mld,'color','r','linestyle','--','LineWidth',1.5)
colortable=textread('MPL_gnuplot.txt');
colormap(f1,colortable);
c=colorbar;
set(c,'ytick',[15 18 21])
caxis([13 23])
ylabel('depth [m]');
text(1.5,-50,'T')
text(17,5,'Transect AB','FontWeight','b')
text(5,-20,'S front','FontWeight','b')
text(21,-35,'transion zone','FontWeight','b')
text(35,-20,'T front','FontWeight','b')
set(gca,'xtick',[],'ytick',[-50 -30 -10])
set(gca,'fontsize',10,'fontweight','b');

f2=axes('Position', [left, bot-zpos*1, width, height]); 
pcolor(x,z,salt);shading interp;colorbar;hold on;
contour(x,z,rho,[22.5:0.15:24],'linewi',.5,'linestyle','-','color','k');
contour(x2,z2,mask1,[1 1],'linewi',1.5,'linestyle','--','color',colorcon);
plot(mld_ABx(1,:),mld,'color','r','linestyle','--','LineWidth',1.5)
colortable=textread('MPL_YlGnBu.txt');
colormap(f2,flipud(colortable));
caxis([30 34.8])
ylabel('depth [m]');
text(1.5,-50,'S')
c=colorbar;
% set(c,'ytick',[31 32.2 33.4])
set(gca,'xtick',[],'ytick',[-50 -30 -10])
set(gca,'fontsize',10,'fontweight','b');


f3=axes('Position', [left, bot-zpos*2, width, height]); 
pcolor(x,z,rho);shading interp;colorbar;hold on;
contour(x,z,rho,[22.5:0.15:24],'linewi',.5,'linestyle','-','color','k');
contour(x2,z2,mask1,[1 1],'linewi',1.5,'linestyle','--','color',colorcon);
plot(mld_ABx(1,:),mld,'color','r','linestyle','--','LineWidth',1.5)
colortable=textread('matlab_jet.txt');
colormap(f3,colortable);
% caxis([22.5 24])
caxis([23 24.5])
c=colorbar;
set(c,'ytick',[23.5 24])
ylabel('depth [m]');
text(1.5,-50,'\sigma')
set(gca,'xtick',[],'ytick',[-50 -30 -10])
set(gca,'fontsize',10,'fontweight','b');


f4=axes('Position', [left, bot-zpos*3, width, height]); 
pcolor(x,z,uacross_AB1);shading interp;colorbar;hold on;
contour(x,z,rho,[22.5:0.15:24],'linewi',.5,'linestyle','-','color','k');
contour(x2,z2,mask1,[1 1],'linewi',1.5,'linestyle','--','color',colorcon);
plot(mld_ABx(1,:),mld,'color','r','linestyle','--','LineWidth',1.5)
colortable=textread('MPL_RdBu.txt');
colormap(f4,flipud(colortable));
caxis([-.5 .5])
ylabel('depth [m]');
text(1.5,-50,'Vel_cross','interpreter','none')
set(gca,'xtick',[],'ytick',[-50 -30 -10])
set(gca,'fontsize',10,'fontweight','b');


f5=axes('Position', [left, bot-zpos*4, width, height]); 
pcolor(x,z,ualong_AB1);shading interp;colorbar;hold on;
contour(x,z,rho,[22.5:0.15:24],'linewi',.5,'linestyle','-','color','k');
contour(x2,z2,mask1,[1 1],'linewi',1.5,'linestyle','--','color',colorcon);
plot(mld_ABx(1,:),mld,'color','r','linestyle','--','LineWidth',1.5)
colortable=textread('MPL_RdBu.txt');
colormap(f5,flipud(colortable));
caxis([-.5 .5])
h = xlabel('km');
set(h, 'Position', [20, -65, -0.5], 'Units', 'normalized');
ylabel('depth [m]');
text(1.5,-50,'Vel_along','interpreter','none')
set(gca,'ytick',[-50 -30 -10],'xtick',[0 10 20 30 40])
set(gca,'fontsize',10,'fontweight','b');
saveas(gcf,'grid100AB','png')



