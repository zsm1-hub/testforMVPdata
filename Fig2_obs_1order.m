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
% [MLDmix,MLDt,MLDr]=get_MLD_obs(temp,rho,z);



%% 原始一阶
figure;
left=0.15;
bot=0.8;
width=0.8;
height=0.17;
zpos=0.18;

% 创建第一个坐标系
f1=axes('Position', [left, bot, width, height]); % 左下角坐标 (0.1, 0.7)，宽度和高度各为 0.25
pcolor(x,z,temp);shading interp;colorbar;hold on;
contour(x,z,rho,[22.5:0.1:24],'linewi',.5,'linestyle','-','color','k');
% plot(x(1,:),MLDt,'color','r')
colortable=textread('MPL_gnuplot.txt');
colormap(f1,colortable);
c=colorbar;
set(c,'ytick',[15 18 21])
caxis([13 23])
ylabel('depth [m]');
text(1.5,-50,'T')
text(17,5,'Transect CD','FontWeight','b')
set(gca,'xtick',[],'ytick',[-50 -30 -10])
set(gca,'fontsize',10,'fontweight','b');

f2=axes('Position', [left, bot-zpos*1, width, height]); 
pcolor(x,z,salt);shading interp;colorbar;hold on;
contour(x,z,rho,[22.5:0.1:24],'linewi',.5,'linestyle','-','color','k');
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

clf
load zsmmvp2.mat
figure;
left=0.15;
bot=0.8;
width=0.8;
height=0.17;
zpos=0.18;

% 创建第一个坐标系
f1=axes('Position', [left, bot, width, height]); % 左下角坐标 (0.1, 0.7)，宽度和高度各为 0.25
pcolor(x,z,temp);shading interp;colorbar;hold on;
contour(x,z,rho,[22.5:0.15:24],'linewi',.5,'linestyle','-','color','k');
% plot(x(1,:),MLDt,'color','r')
colortable=textread('MPL_gnuplot.txt');
colormap(f1,colortable);
c=colorbar;
set(c,'ytick',[15 18 21])
caxis([13 23])
ylabel('depth [m]');
text(1.5,-50,'T')
text(17,5,'Transect AB','FontWeight','b')
set(gca,'xtick',[],'ytick',[-50 -30 -10])
set(gca,'fontsize',10,'fontweight','b');

f2=axes('Position', [left, bot-zpos*1, width, height]); 
pcolor(x,z,salt);shading interp;colorbar;hold on;
contour(x,z,rho,[22.5:0.15:24],'linewi',.5,'linestyle','-','color','k');
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

%% 原始二阶
%% 二阶量

dbdz=v2rho_2d(((rho(1:end-1,:)-rho(2:end,:))./abs(dz))*-g./1025);


%%%thermal wind
uacross2=uacross_CD1;ualong2=ualong_CD1;
dvdz1=dbdx./f;
dvdz=0.5.*(dvdz1(2:end,:)+dvdz1(1:end-1,:));
% dvdz(isnan(dvdz))=0;
%%%%%%%%%%

dvdzreal=v2rho_2d(((uacross2(1:end-1,:)-uacross2(2:end,:))./abs(dz)));
dudzreal=v2rho_2d(((ualong2(1:end-1,:)-ualong2(2:end,:))./abs(dz)));


Ri_balance=(f.^2.*dbdz)./(dbdx.^2);
Ri_gradient=dbdz./(dvdzreal.^2+dudzreal.^2);
dvdx=u2rho_2d((uacross2(:,2:end)-uacross2(:,1:end-1))./dx);
qv=(dvdx+f).*dbdz;qbc=-dvdzreal.*dbdx;q=qv+qbc;
[MLDmix,MLDt,MLDr]=get_MLD_obs(temp,rho,z);
% plot(x(1,:),MLDr,'--','linewi',1.5,'color','b')
% plot(x(1,:),MLDt,'--','linewi',1.5,'color','r')

pycnal=0.2;

f1=subplot(3,3,1)
pcolor(x,z,dvdz1);shading interp;colorbar;hold on;
contour(x,z,rho,[22.5:0.2:24],'linewi',1.2,'linestyle','-','color','k');
plot(x(1,:),MLDmix,'--','linewi',1.5,'color','r')
colortable=textread('NCV_blue_red.txt');
colormap(f1,colortable);
caxis([-.04 .04])
title('dvdz_TW','interpreter','none')
xlabel('km');ylabel('m');
set(gca,'fontsize',12,'fontweight','b');

f2=subplot(3,3,2)
pcolor(x,z,dvdzreal);shading interp;colorbar;hold on;
contour(x,z,rho,[22.5:0.2:24],'linewi',1.2,'linestyle','-','color','k');
plot(x(1,:),MLDr,'--','linewi',1.5,'color','b')
plot(x(1,:),MLDt,'--','linewi',1.5,'color','r')
plot(x(1,:),MLDmix,'--','linewi',1.5,'color','m')

colortable=textread('NCV_blue_red.txt');
colormap(f2,colortable);
caxis([-.04 .04])
title('dvdz','interpreter','none')
xlabel('km');ylabel('m');
set(gca,'fontsize',12,'fontweight','b');

dbdz(dbdz<1e-10)=1e-10;
f3=subplot(3,3,3)
pcolor(x,z,log(dbdz));shading interp;colorbar;hold on;
contour(x,z,rho,[22.5:0.2:24],'linewi',1.2,'linestyle','-','color','k');
plot(x(1,:),MLDr,'--','linewi',1.5,'color','b')
plot(x(1,:),MLDt,'--','linewi',1.5,'color','r')
plot(x(1,:),MLDmix,'--','linewi',1.5,'color','m')

colortable=textread('GMT_ocean.txt');
colormap(f3,flipud(colortable));
caxis([-9 -5])
title('N2','interpreter','none')
xlabel('km');ylabel('m');
set(gca,'fontsize',12,'fontweight','b');

f4=subplot(3,3,4)
pcolor(x,z,Ri_balance);shading interp;colorbar;hold on;
contour(x,z,rho,[22.5:0.2:24],'linewi',1.2,'linestyle','-','color','k');
% contour(x,z,Ri_balance,[0.25 0.25],'linewi',1.2,'linestyle','-','color','w');
% contour(x,z,Ri_balance,[0.75 0.75],'linewi',1.2,'linestyle','-','color','y');
colortable=textread('KH.txt');
colormap(f4,colortable);
c=colorbar;
caxis([0 1.25])
set(c,'ytick',[0.25,1])

title('Ri_balanced','interpreter','none')
xlabel('km');ylabel('m');
set(gca,'fontsize',12,'fontweight','b');

f5=subplot(3,3,5)
pcolor(x,z,Ri_gradient);shading interp;colorbar;hold on;
contour(x,z,rho,[22.5:0.2:24],'linewi',1.2,'linestyle','-','color','k');
% contour(x,z,Ri_gradient,[0.25 0.25],'linewi',1.2,'linestyle','-','color','w');
% contour(x,z,Ri_gradient,[0.75 0.75],'linewi',1.2,'linestyle','-','color','y');
colortable=textread('KH.txt');
colormap(f5,colortable);
c=colorbar;
caxis([0 1.25])
set(c,'ytick',[0.25,1])
title('Ri_balanced','interpreter','none')
xlabel('km');ylabel('m');
set(gca,'fontsize',12,'fontweight','b');

% 
% f6=subplot(3,3,6)
% pcolor(x2./1e3,z2,dtdz2);shading interp;colorbar;hold on;
% contour(x2./1e3,z2,rho2,[22.5:pycnal:24],'linewi',1.2,'linestyle','-','color','k');
% colortable=textread('MPL_RdBu.txt');
% colormap(f6,flipud(colortable));
% caxis([-.2 .2])
% title('dT/dz 500m4m','interpreter','none')
% xlabel('km');ylabel('m');
% set(gca,'fontsize',12,'fontweight','b');

f7=subplot(3,3,7)
pcolor(x,z,qv);shading interp;colorbar;hold on;
contour(x,z,rho,[22.5:0.2:24],'linewi',1.2,'linestyle','-','color','k');
colortable=textread('NCV_blue_red.txt');
colormap(f7,colortable);
caxis([-5e-8 5e-8])
title('qvert','interpreter','none')
xlabel('km');ylabel('m');
set(gca,'fontsize',12,'fontweight','b');

f8=subplot(3,3,8)
pcolor(x,z,qbc);shading interp;colorbar;hold on;
contour(x,z,rho,[22.5:0.2:24],'linewi',1.2,'linestyle','-','color','k');
colortable=textread('NCV_blue_red.txt');
colormap(f8,colortable);
caxis([-5e-8 5e-8])
title('qbc','interpreter','none')
xlabel('km');ylabel('m');
set(gca,'fontsize',12,'fontweight','b');

f9=subplot(3,3,9)
pcolor(x,z,q);shading interp;colorbar;hold on;
contour(x,z,rho,[22.5:0.2:24],'linewi',1.2,'linestyle','-','color','k');
colortable=textread('NCV_blue_red.txt');
colormap(f9,colortable);
caxis([-5e-8 5e-8])
title('q','interpreter','none')
xlabel('km');ylabel('m');
set(gca,'fontsize',12,'fontweight','b');

%%%still too large,because gradient of py is too strong

sgtitle('section CD');
saveas(gcf,'CD2order1','png')

% sgtitle('section AB');
% saveas(gcf,'AB2order1','png')


%% 粗化
xres=500;zres=-2;
xdot=abs(xres)./100;
zdot=abs(zres)./0.5;
clear uacross2;clear ualong2;
for ii=1:floor((size(temp,2)-1)/xdot)
    temp1(:,ii)=nanmean(temp(:,(ii-1)*xdot+1:(ii-1)*xdot+xdot),2);
    salt1(:,ii)=nanmean(salt(:,(ii-1)*xdot+1:(ii-1)*xdot+xdot),2);
    rho1(:,ii)=nanmean(rho(:,(ii-1)*xdot+1:(ii-1)*xdot+xdot),2);
        uacross1(:,ii)=nanmean(uacross_CD1(:,(ii-1)*xdot+1:(ii-1)*xdot+xdot),2);
        ualong1(:,ii)=nanmean(ualong_CD1(:,(ii-1)*xdot+1:(ii-1)*xdot+xdot),2);
    % uacross1(:,ii)=nanmean(uacross_AB1(:,(ii-1)*xdot+1:(ii-1)*xdot+xdot),2);
    % ualong1(:,ii)=nanmean(ualong_AB1(:,(ii-1)*xdot+1:(ii-1)*xdot+xdot),2);
end

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


% [MLDmix,MLDt,MLDr]=get_MLD_obs(temp2,rho2,z2);
% plot(x2(1,:)./1e3,MLDr,'--','linewi',1.5,'color','b')
% plot(x2(1,:)./1e3,MLDt,'--','linewi',1.5,'color','r')

pycnal=0.1;

figure;
left=0.15;
bot=0.8;
width=0.8;
height=0.17;
zpos=0.18;
colorcon='k'

% 创建第一个坐标系
f1=axes('Position', [left, bot, width, height]); % 左下角坐标 (0.1, 0.7)，宽度和高度各为 0.25
pcolor(x2,z2,temp2);shading interp;colorbar;hold on;
contour(x2,z2,rho2,[22.5:pycnal:24],'linewi',.5,'linestyle','-','color',colorcon);
% plot(x(1,:),MLDt,'color','r')
colortable=textread('MPL_gnuplot.txt');
colormap(f1,colortable);
c=colorbar;
set(c,'ytick',[15 18 21])
caxis([13 24])
ylabel('depth [m]');
text(1.5,-50,'T')
text(17,5,'Transect CD','FontWeight','b')
% text(17,5,'Transect AB','FontWeight','b')

set(gca,'xtick',[],'ytick',[-50 -30 -10])
set(gca,'fontsize',10,'fontweight','b');

f2=axes('Position', [left, bot-zpos*1, width, height]); 
pcolor(x2,z2,salt2);shading interp;colorbar;hold on;
contour(x2,z2,rho2,[22.5:pycnal:24],'linewi',.5,'linestyle','-','color',colorcon);
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
colortable=textread('MPL_RdBu.txt');
colormap(f5,flipud(colortable));
caxis([-.5 .5])
h = xlabel('km');
set(h, 'Position', [20, -65, -0.5], 'Units', 'normalized');
ylabel('depth [m]');
text(1.5,-50,'Vel_along','interpreter','none')
set(gca,'ytick',[-50 -30 -10],'xtick',[0 10 20 30 40])
set(gca,'fontsize',10,'fontweight','b');


% saveas(gcf,'grid5001CD','png')
saveas(gcf,'grid5001AB','png')

%% 二阶量

%%%thermal wind

dvdz1=dbdx2./f;
dvdz=0.5.*(dvdz1(2:end,:)+dvdz1(1:end-1,:));dvdz(isnan(dvdz))=0;

dvdzreal=v2rho_2d(((uacross2(1:end-1,:)-uacross2(2:end,:))./abs(zres)));
dudzreal=v2rho_2d(((ualong2(1:end-1,:)-ualong2(2:end,:))./abs(zres)));

effect_rhox_t2=-alpha.*dtdx2;
effect_rhox_s2=beta.*dsdx2;
Rx=-effect_rhox_t2./effect_rhox_s2;
Tu=atan(Rx);

Ri_balance=(f.^2.*dbdz2)./(dbdx2.^2);
Ri_gradient=dbdz2./(dvdzreal.^2+dudzreal.^2);

%%
pycnal=0.1;

figure;
left=0.15;
bot=0.70;
width=0.8;
height=0.25;
zpos=0.3;
colorcon='k'

% 创建第一个坐标系
f1=axes('Position', [left, bot, width, height]); % 左下角坐标 (0.1, 0.7)，宽度和高度各为 0.25
pcolor(x2,z2,dbdx2);shading interp;colorbar;hold on;
contour(x2,z2,rho2,[22.5:pycnal:24],'linewi',.5,'linestyle','-','color',colorcon);
% plot(x(1,:),MLDt,'color','r')
colortable=textread('MPL_RdBu.txt');
colormap(f1,flipud(colortable));
c=colorbar;
% set(c,'ytick',[15 18 21])
caxis([-2e-6 2e-6]);
ylabel('depth [m]');
text(1.5,-50,'M^{2}')
text(17,5,'Transect CD','FontWeight','b')
% text(17,5,'Transect AB','FontWeight','b')
set(gca,'xtick',[],'ytick',[-50 -30 -10])
set(gca,'fontsize',10,'fontweight','b');

f2=axes('Position', [left, bot-zpos*1, width, height]); 
pcolor(x2,z2,log(dbdz2));shading interp;colorbar;hold on;
contour(x2,z2,rho2,[22.5:pycnal:24],'linewi',.5,'linestyle','-','color',colorcon);
colortable=textread('MPL_BuGn.txt');
colormap(f2,(colortable));
caxis([-10 -4])
ylabel('depth [m]');
text(1.5,-50,'N^{2}')
c=colorbar;
% set(c,'ytick',[31 32.2 33.4])
set(gca,'xtick',[],'ytick',[-50 -30 -10])
set(gca,'fontsize',10,'fontweight','b');


f3=axes('Position', [left, bot-zpos*2, width, height]); 
pcolor(x2,z2,Ri_gradient);shading interp;colorbar;hold on;
contour(x2,z2,rho2,[22.5:pycnal:24],'linewi',.5,'linestyle','-','color',colorcon);
colortable=textread('KH.txt');
colormap(f3,colortable);
% caxis([22.5 24])
caxis([0 1.25])
c=colorbar;
set(c,'ytick',[0.25 1])
ylabel('depth [m]');
text(1.5,-50,'Ri_{g}')
set(gca,'ytick',[-50 -30 -10],'xtick',[0 10 20 30 40])
set(gca,'fontsize',10,'fontweight','b');
saveas(gcf,'grid5002CD_2order')
%%


dvdx=u2rho_2d((uacross2(:,2:end)-uacross2(:,1:end-1))./xres);
qv=(dvdx+f).*dbdz2;qbc=-dvdzreal.*dbdx2;q=qv+qbc;

f1=subplot(3,3,1)
pcolor(x2./1e3,z2,dvdz1);shading interp;colorbar;hold on;
contour(x2./1e3,z2,rho2,[22.5:pycnal:24],'linewi',1.2,'linestyle','-','color','k');
% plot(x2(1,:)./1e3,MLDr,'--','linewi',1.5,'color','b')
% plot(x2(1,:)./1e3,MLDt,'--','linewi',1.5,'color','m')
plot(x2(1,:)./1e3,MLDmix,'-','linewi',1.5,'color','b')

colortable=textread('NCV_blue_red.txt');
colormap(f1,colortable);
caxis([-.04 .04])
title('dvdz_TW','interpreter','none')
xlabel('km');ylabel('m');
set(gca,'fontsize',12,'fontweight','b');

f2=subplot(3,3,2)
pcolor(x2./1e3,z2,dvdzreal);shading interp;colorbar;hold on;
contour(x2./1e3,z2,rho2,[22.5:pycnal:24],'linewi',1.2,'linestyle','-','color','k');
% plot(x2(1,:)./1e3,MLDr,'--','linewi',1.5,'color','b')
% plot(x2(1,:)./1e3,MLDt,'--','linewi',1.5,'color','m')
plot(x2(1,:)./1e3,MLDmix,'-','linewi',1.5,'color','b')

colortable=textread('NCV_blue_red.txt');
colormap(f2,colortable);
caxis([-.04 .04])
title('dvdz','interpreter','none')
xlabel('km');ylabel('m');
set(gca,'fontsize',12,'fontweight','b');

dbdz2(dbdz2<1e-10)=1e-10;
f3=subplot(3,3,3)
pcolor(x2./1e3,z2,log(dbdz2));shading interp;colorbar;hold on;
contour(x2./1e3,z2,rho2,[22.5:pycnal:24],'linewi',1.2,'linestyle','-','color','k');
% plot(x2(1,:)./1e3,MLDr,'--','linewi',1.5,'color','b')
% plot(x2(1,:)./1e3,MLDt,'--','linewi',1.5,'color','m')
plot(x2(1,:)./1e3,MLDmix,'-','linewi',1.5,'color','b')

colortable=textread('GMT_ocean.txt');
colormap(f3,flipud(colortable));
caxis([-9 -5])
title('N2','interpreter','none')
xlabel('km');ylabel('m');
set(gca,'fontsize',12,'fontweight','b');

f4=subplot(3,3,4)
pcolor(x2./1e3,z2,Ri_balance);shading interp;colorbar;hold on;
contour(x2./1e3,z2,rho2,[22.5:pycnal:24],'linewi',1.2,'linestyle','-','color','k');
% plot(x2(1,:)./1e3,MLDr,'--','linewi',1.5,'color','b')
% plot(x2(1,:)./1e3,MLDt,'--','linewi',1.5,'color','m')
plot(x2(1,:)./1e3,MLDmix,'-','linewi',1.5,'color','b')

% contour(x2./1e3,z2,Ri_balance,[0.25 0.25],'linewi',1.2,'linestyle','-','color','w');
% contour(x2./1e3,z2,Ri_balance,[0.75 0.75],'linewi',1.2,'linestyle','-','color','y');
colortable=textread('KH.txt');
% colortable=textread('matlab_jet.txt');
colormap(f4,colortable);
c=colorbar;
caxis([0 1.25])
set(c,'ytick',[0.25,1])
caxis([0 10])

title('Ri_balanced','interpreter','none')
xlabel('km');ylabel('m');
set(gca,'fontsize',12,'fontweight','b');

f5=subplot(3,3,5)
pcolor(x2./1e3,z2,Ri_gradient);shading interp;colorbar;hold on;
contour(x2./1e3,z2,rho2,[22.5:pycnal:24],'linewi',1.2,'linestyle','-','color','k');
% plot(x2(1,:)./1e3,MLDr,'--','linewi',1.5,'color','b')
% plot(x2(1,:)./1e3,MLDt,'--','linewi',1.5,'color','m')
plot(x2(1,:)./1e3,MLDmix,'-','linewi',1.5,'color','b')

% contour(x2./1e3,z2,Ri_gradient,[0.25 0.25],'linewi',1.2,'linestyle','-','color','w');
% contour(x2./1e3,z2,Ri_gradient,[0.75 0.75],'linewi',1.2,'linestyle','-','color','y');
colortable=textread('KH.txt');
colormap(f5,colortable);
c=colorbar;
caxis([0 1.25])
set(c,'ytick',[0.25,1])
title('Ri_gradient','interpreter','none')
xlabel('km');ylabel('m');
set(gca,'fontsize',12,'fontweight','b');

% 
% f6=subplot(3,3,6)
% pcolor(x2./1e3,z2,dtdz2);shading interp;colorbar;hold on;
% contour(x2./1e3,z2,rho2,[22.5:pycnal:24],'linewi',1.2,'linestyle','-','color','k');
% colortable=textread('MPL_RdBu.txt');
% colormap(f6,flipud(colortable));
% caxis([-.2 .2])
% title('dT/dz 500m4m','interpreter','none')
% xlabel('km');ylabel('m');
% set(gca,'fontsize',12,'fontweight','b');

f7=subplot(3,3,7)
pcolor(x2./1e3,z2,qv);shading interp;colorbar;hold on;
contour(x2./1e3,z2,rho2,[22.5:pycnal:24],'linewi',1.2,'linestyle','-','color','k');
% plot(x2(1,:)./1e3,MLDr,'--','linewi',1.5,'color','b')
% plot(x2(1,:)./1e3,MLDt,'--','linewi',1.5,'color','m')
plot(x2(1,:)./1e3,MLDmix,'-','linewi',1.5,'color','b')

colortable=textread('NCV_blue_red.txt');
colormap(f7,colortable);
caxis([-5e-8 5e-8])
title('qvert','interpreter','none')
xlabel('km');ylabel('m');
set(gca,'fontsize',12,'fontweight','b');

f8=subplot(3,3,8)
pcolor(x2./1e3,z2,qbc);shading interp;colorbar;hold on;
contour(x2./1e3,z2,rho2,[22.5:pycnal:24],'linewi',1.2,'linestyle','-','color','k');
% plot(x2(1,:)./1e3,MLDr,'--','linewi',1.5,'color','b')
% plot(x2(1,:)./1e3,MLDt,'--','linewi',1.5,'color','m')
plot(x2(1,:)./1e3,MLDmix,'-','linewi',1.5,'color','b')

colortable=textread('NCV_blue_red.txt');
colormap(f8,colortable);
caxis([-5e-8 5e-8])
title('qbc','interpreter','none')
xlabel('km');ylabel('m');
set(gca,'fontsize',12,'fontweight','b');

f9=subplot(3,3,9)
pcolor(x2./1e3,z2,q);shading interp;colorbar;hold on;
contour(x2./1e3,z2,rho2,[22.5:pycnal:24],'linewi',1.2,'linestyle','-','color','k');
% plot(x2(1,:)./1e3,MLDr,'--','linewi',1.5,'color','b')
% plot(x2(1,:)./1e3,MLDt,'--','linewi',1.5,'color','m')
plot(x2(1,:)./1e3,MLDmix,'-','linewi',1.5,'color','b')

colortable=textread('NCV_blue_red.txt');
colormap(f9,colortable);
caxis([-5e-8 5e-8])
title('q','interpreter','none')
xlabel('km');ylabel('m');
set(gca,'fontsize',12,'fontweight','b');

%%%still too large,because gradient of py is too strong

% sgtitle('section CD');
% saveas(gcf,'CD2order1','png')

sgtitle('section AB');
saveas(gcf,'AB2order1','png')

