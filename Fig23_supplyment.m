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
alpha=1.7e-4;
beta=7.6e-4;
cmap = [0 0 1;  % 蓝色 (RGB: 0, 0, 255)
        0 1 0;  % 绿色 (RGB: 0, 255, 0)
        1 1 0;  % 黄色 (RGB: 255, 255, 0)
        1 0 0];

dx=abs(x(1,1)-x(1,2)).*1e3;
dz=abs(z(2,1)-z(1,1));
g=9.81;

dbdz=v2rho_2d(((rho(1:end-1,:)-rho(2:end,:))./abs(dz))*-g./1025);


%%%thermal wind


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

Rx=(alpha.*dtdx2)./(beta.*dsdx2);
Rz=(alpha.*dtdz2)./(beta.*dsdz2);
Tux=atan(Rx);
Tuz=atan(Rz);

% [MLDmix,MLDt,MLDr]=get_MLD_obs(temp2,rho2,z2);
% plot(x2(1,:)./1e3,MLDr,'--','linewi',1.5,'color','b')
% plot(x2(1,:)./1e3,MLDt,'--','linewi',1.5,'color','r')

%%%%MLD PIO
% [mld]=get_mld_taiwanPIO(rho2,z2);
% mld_AB=mld;
% mld_ABx=x2;
% mld_CD=mld;
% mld_CDx=x2;
% save('mld.mat','mld_AB','mld_CD','mld_ABx','mld_CDx');
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
load mld.mat
mld=mld_CD;

% 创建第一个坐标系
f1=axes('Position', [left, bot, width, height]); % 左下角坐标 (0.1, 0.7)，宽度和高度各为 0.25
pcolor(x2,z2,dbdx2);shading interp;colorbar;hold on;
contour(x2,z2,rho2,[22.5:pycnal:24],'linewi',.5,'linestyle','-','color',colorcon,...
    'showtext','on');
plot(mld_CDx(1,:),mld,'color','r','linestyle','--','LineWidth',1.5)
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
mask1=abs(rho2-24.1)<pycnal/1;
%%补充
mask1(1:3,41:44)=1;

pcolor(x2,z2,double(mask1));shading interp;colorbar;hold on;
contour(x2,z2,rho2,[22.5:pycnal:24],'linewi',.5,'linestyle','-','color',colorcon,...
    'showtext','on');
plot(mld_CDx(1,:),mld,'color','r','linestyle','--','LineWidth',1.5)
colortable=textread('MPL_BuGn.txt');
colormap(f2,(colortable));
caxis([0,1])

%%% 分离区域完成
mask_s_front=NaN.*mask1;mask_t_front=NaN.*mask1;
for ii=1:size(z2,1)
    a=mask1(ii,:);
    mask_s_front(ii,1:min(find(a==1)))=1;
    mask_t_front(ii,max(find(a==1)):end)=1;
end
mask_transion=double(mask1);mask_transion(mask_transion==0)=nan;

% N2_S_front=max(max(dbdz2(:,1:38)))
% N2_T_front=max(max(dbdz2(:,39:end)))
% plot(x2(:,1:38),z2(:,1:38))

%% test zone
f3=axes('Position', [left, bot-zpos*2, width, height]); 
pcolor(x2,z2,mask_transion);shading interp;colorbar;hold on;
contour(x2,z2,rho2,[22.5:pycnal:24],'linewi',.5,'linestyle','-','color',colorcon,...
    'showtext','on');
plot(mld_CDx(1,:),mld,'color','r','linestyle','--','LineWidth',1.5)
colortable=textread('MPL_BuGn.txt');
colormap(f3,colortable);
% caxis([22.5 24])
caxis([0 1])
c=colorbar;
set(c,'ytick',[0.25 1])
ylabel('depth [m]');
text(1.5,-50,'Ri_{g}')
set(gca,'ytick',[-50 -30 -10],'xtick',[0 10 20 30 40])
set(gca,'fontsize',10,'fontweight','b');

mask_transion_CD=mask_transion;
mask_t_front_CD=mask_t_front;
mask_s_front_CD=mask_s_front;

max(dbdx2(:))
max(dbdz2(:))

%% 粗化AB
clf
load zsmmvp2.mat;
xdot=abs(xres)./100;
zdot=abs(zres)./0.5;
clear temp1;clear salt1;clear rho1;
clear uacross1;clear ualong1;
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


% [MLDmix,MLDt,MLDr]=get_MLD_obs(temp2,rho2,z2);
% plot(x2(1,:)./1e3,MLDr,'--','linewi',1.5,'color','b')
% plot(x2(1,:)./1e3,MLDt,'--','linewi',1.5,'color','r')

%%%%MLD PIO
% [mld]=get_mld_taiwanPIO(rho2,z2);
% mld_AB=mld;
% mld_ABx=x2;
% save('mld.mat','mld_AB','mld_CD','mld_ABx','mld_CDx');
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
load mld.mat
mld=mld_AB;

% 创建第一个坐标系
f1=axes('Position', [left, bot, width, height]); % 左下角坐标 (0.1, 0.7)，宽度和高度各为 0.25
pcolor(x2,z2,dbdx2);shading interp;colorbar;hold on;
contour(x2,z2,rho2,[22.5:pycnal:25],'linewi',.5,'linestyle','-','color',colorcon,...
    'showtext','on');
plot(mld_ABx(1,:),mld,'color','r','linestyle','--','LineWidth',1.5)
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
mask1=abs(rho2-24.2)<pycnal/1;
%%补充
mask1(1:2,52:63)=1;
pcolor(x2,z2,double(mask1));shading interp;colorbar;hold on;
contour(x2,z2,rho2,[22.5:pycnal:25],'linewi',.5,'linestyle','-','color',colorcon,...
    'showtext','on');
plot(mld_ABx(1,:),mld,'color','r','linestyle','--','LineWidth',1.5)
colortable=textread('MPL_BuGn.txt');
colormap(f2,(colortable));
caxis([0 1])
ylabel('depth [m]');
text(1.5,-50,'N^{2}')
c=colorbar;
% set(c,'ytick',[31 32.2 33.4])
set(gca,'xtick',[],'ytick',[-50 -30 -10])
set(gca,'fontsize',10,'fontweight','b');


%%% 分离区域完成
mask_s_front=NaN.*mask1;mask_t_front=NaN.*mask1;
for ii=1:size(z2,1)
    a=mask1(ii,:);
    mask_s_front(ii,1:min(find(a==1)))=1;
    mask_t_front(ii,max(find(a==1)):end)=1;
end
mask_transion=double(mask1);mask_transion(mask_transion==0)=nan;

% N2_S_front=max(max(dbdz2(:,1:38)))
% N2_T_front=max(max(dbdz2(:,39:end)))
% plot(x2(:,1:38),z2(:,1:38))

%% test zone
f3=axes('Position', [left, bot-zpos*2, width, height]); 
pcolor(x2,z2,mask_s_front);shading interp;colorbar;hold on;
contour(x2,z2,rho2,[22.5:pycnal:25],'linewi',.5,'linestyle','-','color',colorcon,...
    'showtext','on');
plot(mld_ABx(1,:),mld,'color','r','linestyle','--','LineWidth',1.5)
colortable=textread('MPL_BuGn.txt');
colormap(f3,colortable);
% caxis([22.5 24])
caxis([0 1])
c=colorbar;
set(c,'ytick',[0.25 1])
ylabel('depth [m]');
text(1.5,-50,'Ri_{g}')
set(gca,'ytick',[-50 -30 -10],'xtick',[0 10 20 30 40])
set(gca,'fontsize',10,'fontweight','b');

mask_transion_AB=mask_transion;
mask_t_front_AB=mask_t_front;
mask_s_front_AB=mask_s_front;

max(dbdx2(:))
max(dbdz2(:))

save('categories.mat','mask_s_front_AB','mask_t_front_AB','mask_transion_AB',...
    'mask_t_front_CD','mask_s_front_CD','mask_transion_CD');