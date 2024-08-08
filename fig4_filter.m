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

%%%% reality
dvdzreal=v2rho_2d(((uacross2(1:end-1,:)-uacross2(2:end,:))./abs(zres)));
dudzreal=v2rho_2d(((ualong2(1:end-1,:)-ualong2(2:end,:))./abs(zres)));

% effect_rhox_t2=-alpha.*dtdx2;
% effect_rhox_s2=beta.*dsdx2;
% Rx=-effect_rhox_t2./effect_rhox_s2;
% Tu=atan(Rx);

Ri_balance=(f.^2.*dbdz2)./(dbdx2.^2);
Ri_gradient=dbdz2./(dvdzreal.^2+dudzreal.^2);
%% butterworth filter 存在相位差问题
% %% designer a low -pass
% %%% 保留50hz，其他滤掉
% %%频率归一化
% % fs=2*pi*1/500;
% % wp=(2*pi*1/4000)/(fs/2);%4km
% % ws=(2*pi*1/8000)/(fs/2);%8km
% fs=2*pi*1/500;
% wp=2*pi*(1/4000)/(fs/2);%4km
% ws=2*pi*(3.5/4000)/(fs/2);%8km
% alpha_p=3; %%通带允许最大衰减
% alpha_s=60; %%通带允许最小衰减
% %%%%获取阶数和截止频率
% [N1 wc1]=buttord(wp,ws,alpha_p,alpha_s);
% disp(['N1=',num2str(N1),'  wc1=',num2str(wc1)]);
% %%%
% [bb aa]=butter(N1,wc1,'low');
% filter_ucross_CD=NaN.*rho2;
% filter_rho2=NaN.*rho2;
% 
% %%滤波
% for ii=1:16
%     a=~isnan(rho2(ii,:));
%     a1=min(find(a==1));
%     a2=max(find(a==1));
%     if isempty(a1)==1
%         % filter_ucross_CD(ii,a1:a2)=NaN;
%         filter_rho2(ii,a1:a2)=NaN;
% 
%     else
%         % filter_ucross_CD(ii,a1:a2)=filter(bb,aa,uacross2(ii,a1:a2));
%         filter_rho2(ii,a1:a2)=filter(bb,aa,rho2(ii,a1:a2));
%     end
% end
% 
% 
% 
% for ii=1:size(z2,1)
%     b=~isnan(uacross2(ii,:));
%     b1=min(find(b==1));
%     b2=max(find(b==1));
%     if isempty(a1)==1
%         filter_ucross_CD(ii,b1:b2)=NaN;
%     else
%         filter_ucross_CD(ii,b1:b2)=filter(bb,aa,uacross2(ii,b1:b2));
%     end
% end
% 
% % data_interp = interp1(find(~isnan(rho2)), rho2(~isnan(rho2)),...
% %     find(isnan(rho2)), 'linear');
% 
% filter_dbdx2=u2rho_2d((filter_rho2(:,2:end)-filter_rho2(:,1:end-1))./xres)./1025.*-g;
% filter_dvdz1=filter_dbdx2./f;
% 
% filter_dvdzreal=v2rho_2d(((filter_ucross_CD(1:end-1,:)-filter_ucross_CD(2:end,:))./abs(zres)));
%% FIR 滤波
% %%频率归一化
% fs=1/500;
% fc=(1/4000);%4km
% 
% %%%%获取阶数和截止频率
% N=60;
% %%%加窗
% win=hann(N+1);
% b=fir1(N,fc/(fs/2),'low',win);
%% lowpass
fs=2*pi/500;
dxscale=2e3;
% %%滤波
for ii=1:size(z2,1)
    a=~isnan(rho2(ii,:));
    a1=min(find(a==1));
    a2=max(find(a==1));
    if isempty(a1)==1
        % filter_ucross_CD(ii,a1:a2)=NaN;
        filter_rho2(ii,a1:a2)=NaN;

    else
        % filter_ucross_CD(ii,a1:a2)=filter(bb,aa,uacross2(ii,a1:a2));
        % filter_rho2(ii,a1:a2)=filter(bb,aa,rho2(ii,a1:a2));
        filter_rho2(ii,a1:a2)=lowpass(rho2(ii,a1:a2),2*pi/dxscale,fs);

    end
end


uacross2(21:22,1:3)=NaN;
for ii=1:size(z2,1)
    b=~isnan(uacross2(ii,:));
    b1=min(find(b==1));
    b2=max(find(b==1));
    if isempty(b1)==1
        filter_ucross_CD(ii,b1:b2)=NaN;
    else
        % filter_ucross_CD(ii,b1:b2)=filter(bb,aa,uacross2(ii,b1:b2));
        filter_ucross_CD(ii,b1:b2)=lowpass(uacross2(ii,b1:b2),2*pi/dxscale,fs);
    end
end
filter_ucross_CD(filter_ucross_CD==0)=nan;
filter_rho2(filter_rho2==0)=nan;

%%

filter_dbdx2=u2rho_2d((filter_rho2(:,2:end)-filter_rho2(:,1:end-1))./xres)./1025.*-g;
filter_dvdz1=filter_dbdx2./f;

filter_dvdzreal=v2rho_2d(((filter_ucross_CD(1:end-1,:)-filter_ucross_CD(2:end,:))./abs(zres)));
filter_dbdz2=v2rho_2d(((filter_rho2(1:end-1,:)-filter_rho2(2:end,:))./abs(zres))*-g./1025);

filter_Ri_balance=(f.^2.*filter_dbdz2)./(filter_dbdx2.^2);

filter_dvdx=u2rho_2d((filter_ucross_CD(:,2:end)-filter_ucross_CD(:,1:end-1))./xres);
qv=(filter_dvdx+f).*filter_dbdz2;qbc=-filter_dvdzreal.*filter_dbdx2;q=qv+qbc;

%% plot filter CD
pycnal=0.1;

figure;
left=0.15;
bot=0.70;
width=0.4;
height=0.25;
zpos=0.3;
colorcon='k'
load mld.mat
mld=mld_CD;

% mask_CD=mask_transion_CD+mask_s_front_CD.*2+mask_t_front_CD;
mask1=abs(rho2-24.1)<pycnal/1;
%%补充
mask1(1:3,41:44)=1;

% 创建第一个坐标系
f1=axes('Position', [left, bot, width, height]); % 左下角坐标 (0.1, 0.7)，宽度和高度各为 0.25
pcolor(x2,z2,rho2);shading interp;colorbar;hold on;
contour(x2,z2,rho2,[22.5:pycnal:24],'linewi',.5,'linestyle','-','color',colorcon);
% contour(x2,z2,rho2,[24 24],'linewi',1.5,'linestyle','-','color',colorcon);
contour(x2,z2,mask1,[1 1],'linewi',1.5,'linestyle','--','color',colorcon);
plot(mld_CDx(1,:),mld,'color','r','linestyle','--','LineWidth',1.5)
colortable=textread('matlab_jet.txt');
colormap(f1,(colortable));
c=colorbar;
% set(c,'ytick',[15 18 21])
% caxis([-2e-6 2e-6]);
caxis([23 24.5])
ylabel('depth [m]');
text(1.05,-50,'rho')
text(17,5,'Transect CD','FontWeight','b')
set(gca,'xtick',[],'ytick',[-50 -30 -10])
set(gca,'fontsize',10,'fontweight','b');

f2=axes('Position', [left+0.4, bot, width, height]); % 左下角坐标 (0.1, 0.7)，宽度和高度各为 0.25
pcolor(x2,z2,filter_rho2);shading interp;colorbar;hold on;
contour(x2,z2,filter_rho2,[22.5:pycnal:24],'linewi',.5,'linestyle','-','color',colorcon);
% contour(x2,z2,rho2,[24 24],'linewi',1.5,'linestyle','-','color',colorcon);
contour(x2,z2,mask1,[1 1],'linewi',1.5,'linestyle','--','color',colorcon);
plot(mld_CDx(1,:),mld,'color','r','linestyle','--','LineWidth',1.5)
colortable=textread('matlab_jet.txt');
colormap(f2,(colortable));
c=colorbar;
% set(c,'ytick',[15 18 21])
% caxis([-2e-6 2e-6]);
caxis([22.5 23.5])
text(1.05,-50,'f: rho')
set(gca,'xtick',[],'ytick',[])
set(gca,'fontsize',10,'fontweight','b');

f3=axes('Position', [left, bot-zpos*1, width, height]); % 左下角坐标 (0.1, 0.7)，宽度和高度各为 0.25
pcolor(x2,z2,uacross2);shading interp;colorbar;hold on;
contour(x2,z2,rho2,[22.5:pycnal:24],'linewi',.5,'linestyle','-','color',colorcon);
% contour(x2,z2,rho2,[24 24],'linewi',1.5,'linestyle','-','color',colorcon);
contour(x2,z2,mask1,[1 1],'linewi',1.5,'linestyle','--','color',colorcon);
plot(mld_CDx(1,:),mld,'color','r','linestyle','--','LineWidth',1.5)
colortable=textread('MPL_RdBu.txt');
colormap(f3,flipud(colortable));
c=colorbar;
% set(c,'ytick',[15 18 21])
% caxis([-2e-6 2e-6]);
caxis([-.5 .5])
ylabel('depth [m]');
set(gca,'ytick',[-50 -30 -10],'xtick',[0 10 20 30 40])
text(1.05,-50,'ucross')
set(gca,'xtick',[],'ytick',[])
set(gca,'fontsize',10,'fontweight','b');

f4=axes('Position', [left+0.4, bot-zpos*1, width, height]); % 左下角坐标 (0.1, 0.7)，宽度和高度各为 0.25
pcolor(x2,z2,filter_ucross_CD);shading interp;colorbar;hold on;
contour(x2,z2,filter_rho2,[22.5:pycnal:24],'linewi',.5,'linestyle','-','color',colorcon);
% contour(x2,z2,rho2,[24 24],'linewi',1.5,'linestyle','-','color',colorcon);
contour(x2,z2,mask1,[1 1],'linewi',1.5,'linestyle','--','color',colorcon);
plot(mld_CDx(1,:),mld,'color','r','linestyle','--','LineWidth',1.5)
colortable=textread('MPL_RdBu.txt');
colormap(f4,flipud(colortable));
c=colorbar;
% set(c,'ytick',[15 18 21])
% caxis([-2e-6 2e-6]);
caxis([-.5 .5])
text(1.05,-50,'f: ucross')
set(gca,'xtick',[],'ytick',[])
set(gca,'fontsize',10,'fontweight','b');

f5=axes('Position', [left, bot-zpos*2, width, height]); % 左下角坐标 (0.1, 0.7)，宽度和高度各为 0.25
pcolor(x2,z2,filter_dvdz1);shading interp;colorbar;hold on;
contour(x2,z2,filter_rho2,[22.5:pycnal:24],'linewi',.5,'linestyle','-','color',colorcon);
% contour(x2,z2,rho2,[24 24],'linewi',1.5,'linestyle','-','color',colorcon);
contour(x2,z2,mask1,[1 1],'linewi',1.5,'linestyle','--','color',colorcon);
plot(mld_CDx(1,:),mld,'color','r','linestyle','--','LineWidth',1.5)
colortable=textread('MPL_coolwarm.txt');
colormap(f5,(colortable));
c=colorbar;
% set(c,'ytick',[15 18 21])
% caxis([-2e-6 2e-6]);
caxis([-0.025 0.025])
ylabel('depth [m]');
set(gca,'ytick',[],'xtick',[0 10 20 30 40])
ylim([-40 0])
text(1.05,-40,'dv/dz_{TW}')
set(gca,'fontsize',10,'fontweight','b');

f6=axes('Position', [left+0.4, bot-zpos*2, width, height]); % 左下角坐标 (0.1, 0.7)，宽度和高度各为 0.25
pcolor(x2,z2,filter_dvdzreal);shading interp;colorbar;hold on;
contour(x2,z2,filter_rho2,[22.5:pycnal:24],'linewi',.5,'linestyle','-','color',colorcon);
% contour(x2,z2,rho2,[24 24],'linewi',1.5,'linestyle','-','color',colorcon);
contour(x2,z2,mask1,[1 1],'linewi',1.5,'linestyle','--','color',colorcon);
plot(mld_CDx(1,:),mld,'color','r','linestyle','--','LineWidth',1.5)
colortable=textread('MPL_coolwarm.txt');
colormap(f6,(colortable));
c=colorbar;
% set(c,'ytick',[15 18 21])
% caxis([-2e-6 2e-6]);
caxis([-0.025 0.025])
text(1.05,-40,'dv/dz_{data}')
set(gca,'ytick',[],'xtick',[0 10 20 30 40])
ylim([-40 0])
h = xlabel('km');
set(h, 'Position', [-2, -65, -0.5], 'Units', 'normalized');
set(gca,'fontsize',10,'fontweight','b');

saveas(gcf,'fig4_fliter_CD','png')
%% plot filter CD PV
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

% mask_CD=mask_transion_CD+mask_s_front_CD.*2+mask_t_front_CD;
mask1=abs(rho2-24.1)<pycnal/1;
%%补充
mask1(1:3,41:44)=1;

% aaa=filter_Ri_balance(filter_Ri_balance>0);
% filter_Ri_balance(filter_Ri_balance<0)=min(aaa(:));
% % 创建第一个坐标系
% aaa=log(1./filter_Ri_balance);
% aaamask=((aaa<0.5).*(aaa>-0.5));
% % 
% f3=axes('Position', [left, bot-zpos*1, width, height]); % 左下角坐标 (0.1, 0.7)，宽度和高度各为 0.25
% pcolor(x2,z2,log(1./filter_Ri_balance));shading interp;colorbar;hold on;
% % caxis([-1 1])
% contour(x2,z2,rho2,[22.5:pycnal:24],'linewi',.5,'linestyle','-','color',colorcon);
% % contour(x2,z2,rho2,[24 24],'linewi',1.5,'linestyle','-','color',colorcon);
% contour(x2,z2,mask1,[1 1],'linewi',1.5,'linestyle','--','color',colorcon);
% plot(mld_CDx(1,:),mld,'color','r','linestyle','--','LineWidth',1.5)
% colortable=textread('MPL_RdBu.txt');
% colormap(f3,flipud(colortable));
% c=colorbar;
% % set(c,'ytick',[15 18 21])
% % caxis([-2e-6 2e-6]);
% caxis([-1 1])
% ylabel('depth [m]');
% set(gca,'ytick',[-50 -30 -10],'xtick',[0 10 20 30 40])
% text(1.05,-50,'ucross')
% set(gca,'xtick',[],'ytick',[])
% set(gca,'fontsize',10,'fontweight','b');

f4=axes('Position', [left, bot-zpos*0, width, height]); % 左下角坐标 (0.1, 0.7)，宽度和高度各为 0.25
pcolor(x2,z2,qv);shading interp;colorbar;hold on;
contour(x2,z2,filter_rho2,[22.5:pycnal:24],'linewi',.5,'linestyle','-','color',colorcon);
% contour(x2,z2,rho2,[24 24],'linewi',1.5,'linestyle','-','color',colorcon);
contour(x2,z2,mask1,[1 1],'linewi',1.5,'linestyle','--','color',colorcon);
plot(mld_CDx(1,:),mld,'color','r','linestyle','--','LineWidth',1.5)
colortable=textread('MPL_PuOr.txt');
colormap(f4,flipud(colortable));
c=colorbar;
% set(c,'ytick',[15 18 21])
% caxis([-2e-6 2e-6]);
% caxis([-2e-8 2e-8])
caxis([-8e-9 8e-9])

text(17,5,'Transect CD','FontWeight','b')
text(1.05,-50,'f: qv')
set(gca,'xtick',[],'ytick',[-50 -30 -10])
set(gca,'fontsize',10,'fontweight','b');

f5=axes('Position', [left, bot-zpos*1, width, height]); % 左下角坐标 (0.1, 0.7)，宽度和高度各为 0.25
pcolor(x2,z2,qbc);shading interp;colorbar;hold on;
contour(x2,z2,filter_rho2,[22.5:pycnal:24],'linewi',.5,'linestyle','-','color',colorcon);
% contour(x2,z2,rho2,[24 24],'linewi',1.5,'linestyle','-','color',colorcon);
contour(x2,z2,mask1,[1 1],'linewi',1.5,'linestyle','--','color',colorcon);
plot(mld_CDx(1,:),mld,'color','r','linestyle','--','LineWidth',1.5)
colortable=textread('MPL_PuOr.txt');
colormap(f5,flipud(colortable));
c=colorbar;
% set(c,'ytick',[15 18 21])
% caxis([-2e-6 2e-6]);
% caxis([-2e-8 2e-8])
caxis([-8e-9 8e-9])

ylabel('depth [m]');
set(gca,'xtick',[],'ytick',[-50 -30 -10])
text(1.05,-50,'qbc')
set(gca,'fontsize',10,'fontweight','b');

f6=axes('Position', [left, bot-zpos*2, width, height]); % 左下角坐标 (0.1, 0.7)，宽度和高度各为 0.25
pcolor(x2,z2,q);shading interp;colorbar;hold on;
contour(x2,z2,filter_rho2,[22.5:pycnal:24],'linewi',.5,'linestyle','-','color',colorcon);
% contour(x2,z2,rho2,[24 24],'linewi',1.5,'linestyle','-','color',colorcon);
contour(x2,z2,mask1,[1 1],'linewi',1.5,'linestyle','--','color',colorcon);
plot(mld_CDx(1,:),mld,'color','r','linestyle','--','LineWidth',1.5)
colortable=textread('MPL_PuOr.txt');
colormap(f6,flipud(colortable));
c=colorbar;
% set(c,'ytick',[15 18 21])
% caxis([-2e-6 2e-6]);
% caxis([-2e-8 2e-8])
caxis([-8e-9 8e-9])

text(1.05,-50,'q')
set(gca,'xtick',[0 10 20 30 40],'ytick',[-50 -30 -10])
h = xlabel('km');
set(h, 'Position', [20, -65, -0.5], 'Units', 'normalized');
set(gca,'fontsize',10,'fontweight','b');

saveas(gcf,'fig4_fliter_CD_PV','png')
%% instability categories CD
% S_Rib=filter_Ri_balance.*mask_s_front_CD;
% S_xig=filter_dvdx.*mask_s_front_CD;
Rib=filter_Ri_balance;xig=filter_dvdx;N2=filter_dbdz2;
[GIm,GISIm,SIm,ISIm,Stablem]=get_thomas_meassure(N2,Rib,xig,f);

% 
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

% mask_CD=mask_transion_CD+mask_s_front_CD.*2+mask_t_front_CD;
mask1=abs(rho2-24.1)<pycnal/1;
%%补充
mask1(1:3,41:44)=1;
markersize=15;

f1=axes('Position', [left, bot-zpos*0, width, height]); % 左下角坐标 (0.1, 0.7)，宽度和高度各为 0.25
hold on;
scatter(x2.*GIm,z2.*GIm,markersize,[0, 0.5, 0],"filled")
scatter(x2.*GISIm,z2.*GISIm,markersize,[0, 0, 0.5],"filled")
scatter(x2.*SIm,z2.*SIm,markersize,[0.8, 0.6, 0.4],"filled")
scatter(x2.*ISIm,z2.*ISIm,markersize,[0.7, 0.7, 0.7],"filled")
scatter(x2.*Stablem,z2.*Stablem,markersize,'b',"filled")
ax=18;ay=37;
Rib(ax,ay)
xig(ax,ay)
N2(ax,ay)
Rib(ax,ay)-f./(xig(ax,ay)+f)
Stablem(ax,ay)
ISIm(ax,ay)
SIm(ax,ay)
Stablem1=double(N2(ax,ay)>0).*double(xig(ax,ay)<0).*...
    double([Rib(ax,ay)]>(f./(xig(ax,ay)+f))).*double((f./(xig(ax,ay)+f))>1)+...
    double(N2(ax,ay)>0).*double(xig(ax,ay)>0).*double(Rib(ax,ay)>(f./(f+xig(ax,ay))));
double(N2(11,22)>0).*double(xig(11,22)<0).*double(Rib(11,22)>1).*...
    double(Rib<(f./xig(11,22)));
double(N2(11,22)>0).*double(xig(11,22)<0).*...
    double(Rib(11,22)>(f./xig(11,22))).*double((f./xig(11,22))>1)

contour(x2,z2,filter_rho2,[22.5:pycnal:24],'linewi',.5,'linestyle','-','color',colorcon);
contour(x2,z2,mask1,[1 1],'linewi',2,'linestyle','-','color',colorcon);
plot(mld_CDx(1,:),mld,'color','r','linestyle','-','LineWidth',2)
text(17,5,'Transect CD','FontWeight','b')
text(1.05,-50,'f: qv')
set(gca,'xtick',[],'ytick',[-50 -30 -10])
set(gca,'fontsize',10,'fontweight','b');
set(gca,'xtick',[0 10 20 30 40],'ytick',[-50 -30 -10])
h = xlabel('km');
set(h, 'Position', [20, -65, -0.5], 'Units', 'normalized');
set(gca,'fontsize',10,'fontweight','b');

saveas(gcf,'fig4_fliter_CD_Thomas','png')

%% instability and PV CD

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


mask1=abs(rho2-24.1)<pycnal/1;
%%补充
mask1(1:3,41:44)=1;
markersize=15;
f1=axes('Position', [left, bot-zpos*0, width, height]); % 左下角坐标 (0.1, 0.7)，宽度和高度各为 0.25
pcolor(x2,z2,q);shading interp;colorbar;hold on;
contour(x2,z2,filter_rho2,[22.5:pycnal:24],'linewi',.5,'linestyle','-','color',colorcon);
% contour(x2,z2,rho2,[24 24],'linewi',1.5,'linestyle','-','color',colorcon);
contour(x2,z2,mask1,[1 1],'linewi',2,'linestyle','-','color',colorcon);
plot(mld_CDx(1,:),mld,'color','r','linestyle','-','LineWidth',2)
scatter(x2.*GIm,z2.*GIm,markersize,[0, 0.5, 0],"filled")
scatter(x2.*GISIm,z2.*GISIm,markersize,[0, 0, 0.5],"filled")
scatter(x2.*SIm,z2.*SIm,markersize,[0.8, 0.6, 0.4],"filled")
scatter(x2.*ISIm,z2.*ISIm,markersize,[0.7, 0.7, 0.7],"filled")
colortable=textread('MPL_PuOr.txt');
colormap(f1,flipud(colortable));
c=colorbar;


% set(c,'ytick',[15 18 21])
% caxis([-2e-6 2e-6]);
% caxis([-2e-8 2e-8])
caxis([-2e-8 2e-8])
text(1.05,-50,'q')
set(gca,'xtick',[0 10 20 30 40],'ytick',[-50 -30 -10])
h = xlabel('km');
set(h, 'Position', [20, -65, -0.5], 'Units', 'normalized');
set(gca,'fontsize',10,'fontweight','b');

saveas(gcf,'fig4_fliter_CD_THOMAS_PV','png')
%% categories_CD
sky_blue = [135, 206, 235] / 255;
light_green = [144, 238, 144] / 255;
orange = [255, 165, 0] / 255;
mask_current=double(~isnan(q));

S_Rib=filter_Ri_balance.*mask_s_front_CD.*mask_current;
S_xig=(filter_dvdx).*mask_s_front_CD.*mask_current;
S_N2=filter_dbdz2.*mask_s_front_CD.*mask_current;
[S_GIm,S_GISIm,S_SIm,S_ISIm,S_Stablem]=get_thomas_meassure(S_N2,S_Rib,S_xig,f);
total_S=nansum(S_GIm(:))+nansum(S_GISIm(:))+nansum(S_SIm(:))+...
    nansum(S_ISIm(:))+nansum(S_Stablem(:));
S_P_GI=nansum(S_GIm(:))./total_S;
S_P_GISI=nansum(S_GISIm(:))./total_S;
S_P_SI=nansum(S_SIm(:))./total_S;
S_P_ISIm=nansum(S_ISIm(:))./total_S;
S_P_Stable=nansum(S_Stablem(:))./total_S;

S_GIm(isnan(S_GIm))=0;S_GISIm(isnan(S_GISIm))=0;S_SIm(isnan(S_SIm))=0;
S_ISIm(isnan(S_ISIm))=0;S_Stablem(isnan(S_Stablem))=0;

aa=S_GIm+S_GISIm+S_SIm+S_ISIm+S_Stablem
sum(aa(:))
max(aa(:))
[ax,ay]=find(aa==2)




M_Rib=filter_Ri_balance.*mask_transion_CD.*mask_current;
M_xig=filter_dvdx.*mask_transion_CD.*mask_current;
M_N2=filter_dbdz2.*mask_transion_CD.*mask_current;

[M_GIm,M_GISIm,M_SIm,M_ISIm,M_Stablem]=get_thomas_meassure(M_N2,M_Rib,M_xig,f);
total_M=nansum(M_GIm(:))+nansum(M_GISIm(:))+nansum(M_SIm(:))+...
    nansum(M_ISIm(:))+nansum(M_Stablem(:));
M_P_GI=nansum(M_GIm(:))./total_M;
M_P_GISI=nansum(M_GISIm(:))./total_M;
M_P_SI=nansum(M_SIm(:))./total_M;
M_P_ISIm=nansum(M_ISIm(:))./total_M;
M_P_Stable=nansum(M_Stablem(:))./total_M;

% 
% M_GIm(isnan(M_GIm))=0;M_GISIm(isnan(M_GISIm))=0;M_SIm(isnan(M_SIm))=0;
% M_ISIm(isnan(M_ISIm))=0;M_Stablem(isnan(M_Stablem))=0;
% 
% aa=M_GIm+M_GISIm+M_SIm+M_ISIm+M_Stablem
% sum(aa(:))
% max(aa(:))
% [ax,ay]=find(aa==2)



T_Rib=filter_Ri_balance.*mask_t_front_CD.*mask_current;
T_xig=filter_dvdx.*mask_t_front_CD.*mask_current;
T_N2=filter_dbdz2.*mask_t_front_CD.*mask_current;

[T_GIm,T_GISIm,T_SIm,T_ISIm,T_Stablem]=get_thomas_meassure(T_N2,T_Rib,T_xig,f);
total_T=nansum(T_GIm(:))+nansum(T_GISIm(:))+nansum(T_SIm(:))+...
    nansum(T_ISIm(:))+nansum(T_Stablem(:));
T_P_GI=nansum(T_GIm(:))./total_T;
T_P_GISI=nansum(T_GISIm(:))./total_T;
T_P_SI=nansum(T_SIm(:))./total_T;
T_P_ISIm=nansum(T_ISIm(:))./total_T;
T_P_Stable=nansum(T_Stablem(:))./total_T;

% 
% T_GIm(isnan(T_GIm))=0;T_GISIm(isnan(T_GISIm))=0;T_SIm(isnan(T_SIm))=0;
% T_ISIm(isnan(T_ISIm))=0;T_Stablem(isnan(T_Stablem))=0;
% 
% aa=T_GIm+T_GISIm+T_SIm+T_ISIm+T_Stablem
% sum(aa(:))
% max(aa(:))
% [ax,ay]=find(aa==2)


plot([1 2 3 4 5],[S_P_GI,S_P_GISI,S_P_SI,S_P_ISIm,S_P_Stable],'marker','o', ...
    'linestyle','-','LineWidth',2,'color','b');
hold on
plot([1 2 3 4 5],[M_P_GI,M_P_GISI,M_P_SI,M_P_ISIm,M_P_Stable],'marker','o', ...
    'linestyle','-','LineWidth',2,'color','g');
plot([1 2 3 4 5],[T_P_GI,T_P_GISI,T_P_SI,T_P_ISIm,T_P_Stable],'marker','o', ...
    'linestyle','-','LineWidth',2,'color','r');
title('Transect CD','fontsize',10,'fontweight','b')
set(gca,'xtick',[1 2 3 4 5],'XTickLabel',{'GI','GI/SI','SI','CI','Stable'})
ylim([0 0.7])
set(gca,'fontsize',10,'fontweight','b');

saveas(gcf,'fig4_fliter_CD_INST_statisic','png')






%% AB
clear all;close all;clc
load zsmmvp2.mat
load ADCPzsm.mat
load categories.mat

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

%%%% reality
dvdzreal=v2rho_2d(((uacross2(1:end-1,:)-uacross2(2:end,:))./abs(zres)));
dudzreal=v2rho_2d(((ualong2(1:end-1,:)-ualong2(2:end,:))./abs(zres)));

% effect_rhox_t2=-alpha.*dtdx2;
% effect_rhox_s2=beta.*dsdx2;
% Rx=-effect_rhox_t2./effect_rhox_s2;
% Tu=atan(Rx);

Ri_balance=(f.^2.*dbdz2)./(dbdx2.^2);
Ri_gradient=dbdz2./(dvdzreal.^2+dudzreal.^2);
%% butterworth filter 存在相位差问题
% %% designer a low -pass
% %%% 保留50hz，其他滤掉
% %%频率归一化
% % fs=2*pi*1/500;
% % wp=(2*pi*1/4000)/(fs/2);%4km
% % ws=(2*pi*1/8000)/(fs/2);%8km
% fs=2*pi*1/500;
% wp=2*pi*(1/4000)/(fs/2);%4km
% ws=2*pi*(3.5/4000)/(fs/2);%8km
% alpha_p=3; %%通带允许最大衰减
% alpha_s=60; %%通带允许最小衰减
% %%%%获取阶数和截止频率
% [N1 wc1]=buttord(wp,ws,alpha_p,alpha_s);
% disp(['N1=',num2str(N1),'  wc1=',num2str(wc1)]);
% %%%
% [bb aa]=butter(N1,wc1,'low');
% filter_ucross_CD=NaN.*rho2;
% filter_rho2=NaN.*rho2;
% 
% %%滤波
% for ii=1:16
%     a=~isnan(rho2(ii,:));
%     a1=min(find(a==1));
%     a2=max(find(a==1));
%     if isempty(a1)==1
%         % filter_ucross_CD(ii,a1:a2)=NaN;
%         filter_rho2(ii,a1:a2)=NaN;
% 
%     else
%         % filter_ucross_CD(ii,a1:a2)=filter(bb,aa,uacross2(ii,a1:a2));
%         filter_rho2(ii,a1:a2)=filter(bb,aa,rho2(ii,a1:a2));
%     end
% end
% 
% 
% 
% for ii=1:size(z2,1)
%     b=~isnan(uacross2(ii,:));
%     b1=min(find(b==1));
%     b2=max(find(b==1));
%     if isempty(a1)==1
%         filter_ucross_CD(ii,b1:b2)=NaN;
%     else
%         filter_ucross_CD(ii,b1:b2)=filter(bb,aa,uacross2(ii,b1:b2));
%     end
% end
% 
% % data_interp = interp1(find(~isnan(rho2)), rho2(~isnan(rho2)),...
% %     find(isnan(rho2)), 'linear');
% 
% filter_dbdx2=u2rho_2d((filter_rho2(:,2:end)-filter_rho2(:,1:end-1))./xres)./1025.*-g;
% filter_dvdz1=filter_dbdx2./f;
% 
% filter_dvdzreal=v2rho_2d(((filter_ucross_CD(1:end-1,:)-filter_ucross_CD(2:end,:))./abs(zres)));
%% FIR 滤波
% %%频率归一化
% fs=1/500;
% fc=(1/4000);%4km
% 
% %%%%获取阶数和截止频率
% N=60;
% %%%加窗
% win=hann(N+1);
% b=fir1(N,fc/(fs/2),'low',win);
%% lowpass
fs=2*pi/500;
dxscale=4e3;
% %%滤波
for ii=1:size(z2,1)
    a=~isnan(rho2(ii,:));
    a1=min(find(a==1));
    a2=max(find(a==1));
    if isempty(a1)==1
        % filter_ucross_CD(ii,a1:a2)=NaN;
        filter_rho2(ii,a1:a2)=NaN;

    else
        % filter_ucross_CD(ii,a1:a2)=filter(bb,aa,uacross2(ii,a1:a2));
        % filter_rho2(ii,a1:a2)=filter(bb,aa,rho2(ii,a1:a2));
        filter_rho2(ii,a1:a2)=lowpass(rho2(ii,a1:a2),2*pi/dxscale,fs);

    end
end


uacross2(21:22,1:3)=NaN;
for ii=1:size(z2,1)
    b=~isnan(uacross2(ii,:));
    b1=min(find(b==1));
    b2=max(find(b==1));
    if isempty(b1)==1
        filter_ucross_AB(ii,b1:b2)=NaN;
    else
        % filter_ucross_CD(ii,b1:b2)=filter(bb,aa,uacross2(ii,b1:b2));
        filter_ucross_AB(ii,b1:b2)=lowpass(uacross2(ii,b1:b2),2*pi/dxscale,fs);
    end
end
filter_ucross_AB(filter_ucross_AB==0)=nan;
filter_rho2(filter_rho2==0)=nan;

%%

filter_dbdx2=u2rho_2d((filter_rho2(:,2:end)-filter_rho2(:,1:end-1))./xres)./1025.*-g;
filter_dvdz1=filter_dbdx2./f;

filter_dvdzreal=v2rho_2d(((filter_ucross_AB(1:end-1,:)-filter_ucross_AB(2:end,:))./abs(zres)));
filter_dbdz2=v2rho_2d(((filter_rho2(1:end-1,:)-filter_rho2(2:end,:))./abs(zres))*-g./1025);

filter_Ri_balance=(f.^2.*filter_dbdz2)./(filter_dbdx2.^2);

filter_dvdx=u2rho_2d((filter_ucross_AB(:,2:end)-filter_ucross_AB(:,1:end-1))./xres);
qv=(filter_dvdx+f).*filter_dbdz2;qbc=-filter_dvdzreal.*filter_dbdx2;q=qv+qbc;

%% 
pycnal=0.1;

figure;
left=0.15;
bot=0.70;
width=0.4;
height=0.25;
zpos=0.3;
colorcon='k'
load mld.mat
mld=mld_AB;

mask1=abs(rho2-24.2)<pycnal/1;
%%补充
mask1(1:2,52:63)=1;

% 创建第一个坐标系
f1=axes('Position', [left, bot, width, height]); % 左下角坐标 (0.1, 0.7)，宽度和高度各为 0.25
pcolor(x2,z2,rho2);shading interp;colorbar;hold on;
contour(x2,z2,rho2,[22.5:pycnal:24],'linewi',.5,'linestyle','-','color',colorcon);
% contour(x2,z2,rho2,[24 24],'linewi',1.5,'linestyle','-','color',colorcon);
contour(x2,z2,mask1,[1 1],'linewi',1.5,'linestyle','--','color',colorcon);
plot(mld_ABx(1,:),mld,'color','r','linestyle','--','LineWidth',1.5)
colortable=textread('matlab_jet.txt');
colormap(f1,(colortable));
c=colorbar;
% set(c,'ytick',[15 18 21])
% caxis([-2e-6 2e-6]);
caxis([23 24.5])
ylabel('depth [m]');
text(1.05,-50,'rho')
text(17,5,'Transect AB','FontWeight','b')
set(gca,'xtick',[],'ytick',[-50 -30 -10])
set(gca,'fontsize',10,'fontweight','b');

f2=axes('Position', [left+0.4, bot, width, height]); % 左下角坐标 (0.1, 0.7)，宽度和高度各为 0.25
pcolor(x2,z2,filter_rho2);shading interp;colorbar;hold on;
contour(x2,z2,filter_rho2,[22.5:pycnal:24],'linewi',.5,'linestyle','-','color',colorcon);
% contour(x2,z2,rho2,[24 24],'linewi',1.5,'linestyle','-','color',colorcon);
contour(x2,z2,mask1,[1 1],'linewi',1.5,'linestyle','--','color',colorcon);
plot(mld_ABx(1,:),mld,'color','r','linestyle','--','LineWidth',1.5)
colortable=textread('matlab_jet.txt');
colormap(f2,(colortable));
c=colorbar;
% set(c,'ytick',[15 18 21])
% caxis([-2e-6 2e-6]);
caxis([22.5 23.5])
text(1.05,-50,'f: rho')
set(gca,'xtick',[],'ytick',[])
set(gca,'fontsize',10,'fontweight','b');

f3=axes('Position', [left, bot-zpos*1, width, height]); % 左下角坐标 (0.1, 0.7)，宽度和高度各为 0.25
pcolor(x2,z2,uacross2);shading interp;colorbar;hold on;
contour(x2,z2,rho2,[22.5:pycnal:24],'linewi',.5,'linestyle','-','color',colorcon);
% contour(x2,z2,rho2,[24 24],'linewi',1.5,'linestyle','-','color',colorcon);
contour(x2,z2,mask1,[1 1],'linewi',1.5,'linestyle','--','color',colorcon);
plot(mld_ABx(1,:),mld,'color','r','linestyle','--','LineWidth',1.5)
colortable=textread('MPL_RdBu.txt');
colormap(f3,flipud(colortable));
c=colorbar;
% set(c,'ytick',[15 18 21])
% caxis([-2e-6 2e-6]);
caxis([-.5 .5])
ylabel('depth [m]');
set(gca,'ytick',[-50 -30 -10],'xtick',[0 10 20 30 40])
text(1.05,-50,'ucross')
set(gca,'xtick',[],'ytick',[])
set(gca,'fontsize',10,'fontweight','b');

f4=axes('Position', [left+0.4, bot-zpos*1, width, height]); % 左下角坐标 (0.1, 0.7)，宽度和高度各为 0.25
pcolor(x2,z2,filter_ucross_AB);shading interp;colorbar;hold on;
contour(x2,z2,filter_rho2,[22.5:pycnal:24],'linewi',.5,'linestyle','-','color',colorcon);
% contour(x2,z2,rho2,[24 24],'linewi',1.5,'linestyle','-','color',colorcon);
contour(x2,z2,mask1,[1 1],'linewi',1.5,'linestyle','--','color',colorcon);
plot(mld_ABx(1,:),mld,'color','r','linestyle','--','LineWidth',1.5)
colortable=textread('MPL_RdBu.txt');
colormap(f4,flipud(colortable));
c=colorbar;
% set(c,'ytick',[15 18 21])
% caxis([-2e-6 2e-6]);
caxis([-.5 .5])
text(1.05,-50,'f: ucross')
set(gca,'xtick',[],'ytick',[])
set(gca,'fontsize',10,'fontweight','b');

f5=axes('Position', [left, bot-zpos*2, width, height]); % 左下角坐标 (0.1, 0.7)，宽度和高度各为 0.25
pcolor(x2,z2,filter_dvdz1);shading interp;colorbar;hold on;
contour(x2,z2,filter_rho2,[22.5:pycnal:24],'linewi',.5,'linestyle','-','color',colorcon);
% contour(x2,z2,rho2,[24 24],'linewi',1.5,'linestyle','-','color',colorcon);
contour(x2,z2,mask1,[1 1],'linewi',1.5,'linestyle','--','color',colorcon);
plot(mld_ABx(1,:),mld,'color','r','linestyle','--','LineWidth',1.5)
colortable=textread('MPL_coolwarm.txt');
colormap(f5,(colortable));
c=colorbar;
% set(c,'ytick',[15 18 21])
% caxis([-2e-6 2e-6]);
caxis([-0.025 0.025])
ylabel('depth [m]');
set(gca,'ytick',[-50 -30 -10],'xtick',[0 10 20 30 40])
text(1.05,-40,'dv/dz_{TW}')
ylim([-40 0])
set(gca,'fontsize',10,'fontweight','b');

f6=axes('Position', [left+0.4, bot-zpos*2, width, height]); % 左下角坐标 (0.1, 0.7)，宽度和高度各为 0.25
pcolor(x2,z2,filter_dvdzreal);shading interp;colorbar;hold on;
contour(x2,z2,filter_rho2,[22.5:pycnal:24],'linewi',.5,'linestyle','-','color',colorcon);
% contour(x2,z2,rho2,[24 24],'linewi',1.5,'linestyle','-','color',colorcon);
contour(x2,z2,mask1,[1 1],'linewi',1.5,'linestyle','--','color',colorcon);
plot(mld_ABx(1,:),mld,'color','r','linestyle','--','LineWidth',1.5)
colortable=textread('MPL_coolwarm.txt');
colormap(f6,(colortable));
c=colorbar;
% set(c,'ytick',[15 18 21])
% caxis([-2e-6 2e-6]);
caxis([-0.025 0.025])
text(1.05,-40,'dv/dz_{data}')
set(gca,'ytick',[],'xtick',[0 10 20 30 40])
ylim([-40 0])
h = xlabel('km');
set(h, 'Position', [-2, -65, -0.5], 'Units', 'normalized');
set(gca,'fontsize',10,'fontweight','b');


saveas(gcf,'fig4_fliter_AB','png')

%% plot filter AB PV
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


mask1=abs(rho2-24.2)<pycnal/1;
%%补充
mask1(1:2,52:63)=1;

% aaa=filter_Ri_balance(filter_Ri_balance>0);
% filter_Ri_balance(filter_Ri_balance<0)=min(aaa(:));
% % 创建第一个坐标系
% aaa=log(1./filter_Ri_balance);
% aaamask=((aaa<0.5).*(aaa>-0.5));
% % 
% f3=axes('Position', [left, bot-zpos*1, width, height]); % 左下角坐标 (0.1, 0.7)，宽度和高度各为 0.25
% pcolor(x2,z2,log(1./filter_Ri_balance));shading interp;colorbar;hold on;
% % caxis([-1 1])
% contour(x2,z2,rho2,[22.5:pycnal:24],'linewi',.5,'linestyle','-','color',colorcon);
% % contour(x2,z2,rho2,[24 24],'linewi',1.5,'linestyle','-','color',colorcon);
% contour(x2,z2,mask1,[1 1],'linewi',1.5,'linestyle','--','color',colorcon);
% plot(mld_CDx(1,:),mld,'color','r','linestyle','--','LineWidth',1.5)
% colortable=textread('MPL_RdBu.txt');
% colormap(f3,flipud(colortable));
% c=colorbar;
% % set(c,'ytick',[15 18 21])
% % caxis([-2e-6 2e-6]);
% caxis([-1 1])
% ylabel('depth [m]');
% set(gca,'ytick',[-50 -30 -10],'xtick',[0 10 20 30 40])
% text(1.05,-50,'ucross')
% set(gca,'xtick',[],'ytick',[])
% set(gca,'fontsize',10,'fontweight','b');

f4=axes('Position', [left, bot-zpos*0, width, height]); % 左下角坐标 (0.1, 0.7)，宽度和高度各为 0.25
pcolor(x2,z2,qv);shading interp;colorbar;hold on;
contour(x2,z2,filter_rho2,[22.5:pycnal:24],'linewi',.5,'linestyle','-','color',colorcon);
% contour(x2,z2,rho2,[24 24],'linewi',1.5,'linestyle','-','color',colorcon);
contour(x2,z2,mask1,[1 1],'linewi',1.5,'linestyle','--','color',colorcon);
plot(mld_ABx(1,:),mld,'color','r','linestyle','--','LineWidth',1.5)
colortable=textread('MPL_PuOr.txt');
colormap(f4,flipud(colortable));
c=colorbar;
% set(c,'ytick',[15 18 21])
% caxis([-2e-6 2e-6]);
% caxis([-2e-8 2e-8])
caxis([-2e-9 2e-9])

text(17,5,'Transect AB','FontWeight','b')
text(1.05,-50,'f: qv')
set(gca,'xtick',[],'ytick',[-50 -30 -10])
set(gca,'fontsize',10,'fontweight','b');

f5=axes('Position', [left, bot-zpos*1, width, height]); % 左下角坐标 (0.1, 0.7)，宽度和高度各为 0.25
pcolor(x2,z2,qbc);shading interp;colorbar;hold on;
contour(x2,z2,filter_rho2,[22.5:pycnal:24],'linewi',.5,'linestyle','-','color',colorcon);
% contour(x2,z2,rho2,[24 24],'linewi',1.5,'linestyle','-','color',colorcon);
contour(x2,z2,mask1,[1 1],'linewi',1.5,'linestyle','--','color',colorcon);
plot(mld_ABx(1,:),mld,'color','r','linestyle','--','LineWidth',1.5)
colortable=textread('MPL_PuOr.txt');
colormap(f5,flipud(colortable));
c=colorbar;
% set(c,'ytick',[15 18 21])
% caxis([-2e-6 2e-6]);
% caxis([-2e-8 2e-8])
caxis([-2e-9 2e-9])

ylabel('depth [m]');
set(gca,'xtick',[],'ytick',[-50 -30 -10])
text(1.05,-50,'qbc')
set(gca,'fontsize',10,'fontweight','b');

f6=axes('Position', [left, bot-zpos*2, width, height]); % 左下角坐标 (0.1, 0.7)，宽度和高度各为 0.25
pcolor(x2,z2,q);shading interp;colorbar;hold on;
contour(x2,z2,filter_rho2,[22.5:pycnal:24],'linewi',.5,'linestyle','-','color',colorcon);
% contour(x2,z2,rho2,[24 24],'linewi',1.5,'linestyle','-','color',colorcon);
contour(x2,z2,mask1,[1 1],'linewi',1.5,'linestyle','--','color',colorcon);
plot(mld_ABx(1,:),mld,'color','r','linestyle','--','LineWidth',1.5)
colortable=textread('MPL_PuOr.txt');
colormap(f6,flipud(colortable));
c=colorbar;
% set(c,'ytick',[15 18 21])
% caxis([-2e-6 2e-6]);
% caxis([-2e-8 2e-8])
caxis([-2e-9 2e-9])
text(1.05,-50,'q')
set(gca,'xtick',[0 10 20 30 40],'ytick',[-50 -30 -10])
h = xlabel('km');
set(h, 'Position', [20, -65, -0.5], 'Units', 'normalized');
set(gca,'fontsize',10,'fontweight','b');

saveas(gcf,'fig4_fliter_AB_PV','png')


%% instability categories AB
% S_Rib=filter_Ri_balance.*mask_s_front_CD;
% S_xig=filter_dvdx.*mask_s_front_CD;
Rib=filter_Ri_balance;xig=filter_dvdx;N2=filter_dbdz2;

[GIm,GISIm,SIm,ISIm,Stablem]=get_thomas_meassure(N2,Rib,xig,f);

% 
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

mask1=abs(rho2-24.2)<pycnal/1;
%%补充
mask1(1:2,52:63)=1;
markersize=15;

f1=axes('Position', [left, bot-zpos*0, width, height]); % 左下角坐标 (0.1, 0.7)，宽度和高度各为 0.25
hold on;
scatter(x2.*GIm,z2.*GIm,markersize,[0, 0.5, 0],"filled")
scatter(x2.*GISIm,z2.*GISIm,markersize,[0, 0, 0.5],"filled")
scatter(x2.*SIm,z2.*SIm,markersize,[0.8, 0.6, 0.4],"filled")
scatter(x2.*ISIm,z2.*ISIm,markersize,[0.7, 0.7, 0.7],"filled")
scatter(x2.*Stablem,z2.*Stablem,markersize,'b',"filled")

contour(x2,z2,filter_rho2,[22.5:pycnal:24],'linewi',.5,'linestyle','-','color',colorcon);
contour(x2,z2,mask1,[1 1],'linewi',2,'linestyle','-','color',colorcon);
plot(mld_ABx(1,:),mld,'color','r','linestyle','-','LineWidth',2)
text(17,5,'Transect AB','FontWeight','b')
text(1.05,-50,'f: qv')
set(gca,'xtick',[],'ytick',[-50 -30 -10])
set(gca,'fontsize',10,'fontweight','b');
set(gca,'xtick',[0 10 20 30 40],'ytick',[-50 -30 -10])
h = xlabel('km');
set(h, 'Position', [20, -65, -0.5], 'Units', 'normalized');
set(gca,'fontsize',10,'fontweight','b');

saveas(gcf,'fig4_fliter_AB_Thomas','png')
%% instability and PV

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


mask1=abs(rho2-24.2)<pycnal/1;
%%补充
mask1(1:2,52:63)=1;
f1=axes('Position', [left, bot-zpos*0, width, height]); % 左下角坐标 (0.1, 0.7)，宽度和高度各为 0.25
pcolor(x2,z2,q);shading interp;colorbar;hold on;
contour(x2,z2,filter_rho2,[22.5:pycnal:24],'linewi',.5,'linestyle','-','color',colorcon);
% contour(x2,z2,rho2,[24 24],'linewi',1.5,'linestyle','-','color',colorcon);
contour(x2,z2,mask1,[1 1],'linewi',2,'linestyle','-','color',colorcon);
plot(mld_ABx(1,:),mld,'color','r','linestyle','-','LineWidth',2)
scatter(x2.*GIm,z2.*GIm,markersize,[0, 0.5, 0],"filled")
scatter(x2.*GISIm,z2.*GISIm,markersize,[0, 0, 0.5],"filled")
scatter(x2.*SIm,z2.*SIm,markersize,[0.8, 0.6, 0.4],"filled")
scatter(x2.*ISIm,z2.*ISIm,markersize,[0.7, 0.7, 0.7],"filled")
colortable=textread('MPL_PuOr.txt');
colormap(f1,flipud(colortable));
c=colorbar;


% set(c,'ytick',[15 18 21])
% caxis([-2e-6 2e-6]);
% caxis([-2e-8 2e-8])
caxis([-2e-9 2e-9])
text(1.05,-50,'q')
set(gca,'xtick',[0 10 20 30 40],'ytick',[-50 -30 -10])
h = xlabel('km');
set(h, 'Position', [20, -65, -0.5], 'Units', 'normalized');
set(gca,'fontsize',10,'fontweight','b');

saveas(gcf,'fig4_fliter_AB_THOMAS_PV','png')



%% categories_AB
sky_blue = [135, 206, 235] / 255;
light_green = [144, 238, 144] / 255;
orange = [255, 165, 0] / 255;
mask_current=double(~isnan(q));

S_Rib=filter_Ri_balance.*mask_s_front_AB.*mask_current;
S_xig=filter_dvdx.*mask_s_front_AB.*mask_current;
S_N2=filter_dbdz2.*mask_s_front_AB.*mask_current;

[S_GIm,S_GISIm,S_SIm,S_ISIm,S_Stablem]=get_thomas_meassure(S_N2,S_Rib,S_xig,f);
total_S=nansum(S_GIm(:))+nansum(S_GISIm(:))+nansum(S_SIm(:))+...
    nansum(S_ISIm(:))+nansum(S_Stablem(:));
S_P_GI=nansum(S_GIm(:))./total_S;
S_P_GISI=nansum(S_GISIm(:))./total_S;
S_P_SI=nansum(S_SIm(:))./total_S;
S_P_ISIm=nansum(S_ISIm(:))./total_S;
S_P_Stable=nansum(S_Stablem(:))./total_S;

% S_GIm(isnan(S_GIm))=0;S_GISIm(isnan(S_GISIm))=0;S_SIm(isnan(S_SIm))=0;
% S_ISIm(isnan(S_ISIm))=0;S_Stablem(isnan(S_Stablem))=0;
% 
% aa=S_GIm+S_GISIm+S_SIm+S_ISIm+S_Stablem
% sum(aa(:))
% max(aa(:))
% [ax,ay]=find(aa==2)


M_Rib=filter_Ri_balance.*mask_transion_AB.*mask_current;
M_xig=filter_dvdx.*mask_transion_AB.*mask_current;
M_N2=filter_dbdz2.*mask_transion_AB.*mask_current;

[M_GIm,M_GISIm,M_SIm,M_ISIm,M_Stablem]=get_thomas_meassure(M_N2,M_Rib,M_xig,f);
total_M=nansum(M_GIm(:))+nansum(M_GISIm(:))+nansum(M_SIm(:))+...
    nansum(M_ISIm(:))+nansum(M_Stablem(:));
M_P_GI=nansum(M_GIm(:))./total_M;
M_P_GISI=nansum(M_GISIm(:))./total_M;
M_P_SI=nansum(M_SIm(:))./total_M;
M_P_ISIm=nansum(M_ISIm(:))./total_M;
M_P_Stable=nansum(M_Stablem(:))./total_M;

% M_GIm(isnan(M_GIm))=0;M_GISIm(isnan(M_GISIm))=0;M_SIm(isnan(M_SIm))=0;
% M_ISIm(isnan(M_ISIm))=0;M_Stablem(isnan(M_Stablem))=0;
% 
% aa=M_GIm+M_GISIm+M_SIm+M_ISIm+M_Stablem
% sum(aa(:))
% max(aa(:))
% [ax,ay]=find(aa==2)

T_Rib=filter_Ri_balance.*mask_t_front_AB.*mask_current;
T_xig=filter_dvdx.*mask_t_front_AB.*mask_current;
T_N2=filter_dbdz2.*mask_t_front_AB.*mask_current;

[T_GIm,T_GISIm,T_SIm,T_ISIm,T_Stablem]=get_thomas_meassure(T_N2,T_Rib,T_xig,f);
total_T=nansum(T_GIm(:))+nansum(T_GISIm(:))+nansum(T_SIm(:))+...
    nansum(T_ISIm(:))+nansum(T_Stablem(:));
T_P_GI=nansum(T_GIm(:))./total_T;
T_P_GISI=nansum(T_GISIm(:))./total_T;
T_P_SI=nansum(T_SIm(:))./total_T;
T_P_ISIm=nansum(T_ISIm(:))./total_T;
T_P_Stable=nansum(T_Stablem(:))./total_T;



% T_GIm(isnan(T_GIm))=0;T_GISIm(isnan(T_GISIm))=0;T_SIm(isnan(T_SIm))=0;
% T_ISIm(isnan(T_ISIm))=0;T_Stablem(isnan(T_Stablem))=0;
% 
% aa=T_GIm+T_GISIm+T_SIm+T_ISIm+T_Stablem
% sum(aa(:))
% max(aa(:))
% [ax,ay]=find(aa==2)


plot([1 2 3 4 5],[S_P_GI,S_P_GISI,S_P_SI,S_P_ISIm,S_P_Stable],'marker','o', ...
    'linestyle','-','LineWidth',2,'color',sky_blue);
hold on
plot([1 2 3 4 5],[M_P_GI,M_P_GISI,M_P_SI,M_P_ISIm,M_P_Stable],'marker','o', ...
    'linestyle','-','LineWidth',2,'color',light_green);
plot([1 2 3 4 5],[T_P_GI,T_P_GISI,T_P_SI,T_P_ISIm,T_P_Stable],'marker','o', ...
    'linestyle','-','LineWidth',2,'color',orange);
title('Transect AB','fontsize',10,'fontweight','b')
set(gca,'xtick',[1 2 3 4 5],'XTickLabel',{'GI','GI/SI','SI','CI','Stable'})
set(gca,'fontsize',10,'fontweight','b');
ylim([0 0.7])
saveas(gcf,'fig4_fliter_AB_INST_statisic','png')
% saveas(gcf,'fig4_fliter_ALL_INST_statisic','png')
