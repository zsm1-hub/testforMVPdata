clear all;close all;clc
addpath(genpath('E:\ROMS学习\download_data_process\submeso\analysis\GSW\seawater\seawater'));
addpath('F:\TWS_Acrobat\TWS_Acrobat\TWS_Acrobat\')
addpath('E:\ROMS学习\download_data_process\submeso\initial')
addpath('E:\ROMS学习\download_data_process\submeso\analysis\taiwan')
% addpath('E:\ROMS学习\download_data_process\colorbar\colorbar_NCL\colorbar_NCL')
addpath('D:\colorbar\colorbar_NCL');
% load zsmmvp1.mat
load zsmmvp2.mat
load ADCPzsm.mat
%%%%1是CD/FG，2是AB/DE

dx=abs(x(1,1)-x(1,2)).*1e3;
dz=abs(z(2,1)-z(1,1));
g=9.81;

dbdz=v2rho_2d(((rho(1:end-1,:)-rho(2:end,:))./abs(dz))*-g./1025);


%%%thermal wind
uacross2=uacross_CD1;ualong2=ualong_CD1;
% dvdz1=dbdx./f;
% dvdz=0.5.*(dvdz1(2:end,:)+dvdz1(1:end-1,:));
% dvdz(isnan(dvdz))=0;
%%%%%%%%%%

% dvdzreal=v2rho_2d(((uacross2(1:end-1,:)-uacross2(2:end,:))./abs(dz)));
% dudzreal=v2rho_2d(((ualong2(1:end-1,:)-ualong2(2:end,:))./abs(dz)));
% 
% 
% Ri_balance=(f.^2.*dbdz)./(dbdx.^2);
% Ri_gradient=dbdz./(dvdzreal.^2+dudzreal.^2);

% dvdx=u2rho_2d((uacross2(:,2:end)-uacross2(:,1:end-1))./dx);
% qv=(dvdx+f).*dbdz;qbc=-dvdzreal.*dbdx;q=qv+qbc;
% [MLDmix,MLDt,MLDr]=get_MLD_obs(temp,rho,z);
% plot(x(1,:),MLDr,'--','linewi',1.5,'color','b')
% plot(x(1,:),MLDt,'--','linewi',1.5,'color','r')

%% 粗化
xres=500;zres=-2;
xdot=abs(xres)./100;
zdot=abs(zres)./0.5;
clear uacross2;clear ualong2;
for ii=1:floor((size(temp,2)-1)/xdot)
    temp1(:,ii)=nanmean(temp(:,(ii-1)*xdot+1:(ii-1)*xdot+xdot),2);
    salt1(:,ii)=nanmean(salt(:,(ii-1)*xdot+1:(ii-1)*xdot+xdot),2);
    rho1(:,ii)=nanmean(rho(:,(ii-1)*xdot+1:(ii-1)*xdot+xdot),2);
        % uacross1(:,ii)=nanmean(uacross_CD1(:,(ii-1)*xdot+1:(ii-1)*xdot+xdot),2);
        % ualong1(:,ii)=nanmean(ualong_CD1(:,(ii-1)*xdot+1:(ii-1)*xdot+xdot),2);
    uacross1(:,ii)=nanmean(uacross_AB1(:,(ii-1)*xdot+1:(ii-1)*xdot+xdot),2);
    ualong1(:,ii)=nanmean(ualong_AB1(:,(ii-1)*xdot+1:(ii-1)*xdot+xdot),2);
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



%% 二阶量

%%%thermal wind

dvdz1=dbdx2./f;
dvdz=0.5.*(dvdz1(2:end,:)+dvdz1(1:end-1,:));dvdz(isnan(dvdz))=0;

dvdzreal=v2rho_2d(((uacross2(1:end-1,:)-uacross2(2:end,:))./abs(zres)));
dudzreal=v2rho_2d(((ualong2(1:end-1,:)-ualong2(2:end,:))./abs(zres)));

% effect_rhox_t2=-alpha.*dtdx2;
% effect_rhox_s2=beta.*dsdx2;
% Rx=-effect_rhox_t2./effect_rhox_s2;
% Tu=atan(Rx);

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
% text(17,5,'Transect CD','FontWeight','b')
text(17,5,'Transect AB','FontWeight','b')
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

% saveas(gcf,'grid5002CD_2order','png')
saveas(gcf,'grid5002AB_2order','png')


max(dbdx2(:))
max(dbdz2(:))