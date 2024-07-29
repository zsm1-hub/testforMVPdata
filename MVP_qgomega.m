clear all;close all;clc

f=gsw_f(24.4);g=9.8;
load mvppos.mat
load mvpadcp.mat
dx=abs(Ds_low(1,2)-Ds_low(1,1)).*1e3;
dz=abs(De_low(2,1)-De_low(1,1));
%%%%%假设混合层内地转流均一
ualong2=repmat(ualong1,[size(Ds_low,1) 1]);
rho1=PotDen_low+1000;
b=g.*(1-rho1./1025);
dbdz1=(b(2:end,:)-b(1:end-1,:))./dz;
N2=v2rho_2d(dbdz1);
Nthreshold=min(abs(N2(:)));
N2(N2==0)=Nthreshold;
N2(N2<0)=Nthreshold;


pcolor(Ds_low,De_low,N2);shading flat;colorbar;hold on
% colortable = textread('GMT_polar.txt');
colortable = textread('MPL_viridis.txt');
contour(Ds_low,De_low,PotDen_low,'color','w','linewi',1.4);
colormap(colortable)
caxis([0 5e-4])
title('N2');
%%

clear all;close all;clc

f=gsw_f(24.4);g=9.8;
load mvppos.mat
load mvpadcp.mat
%for ED
% xrange=40:431;
% yrange=13:82;
% for FG
xrange=125:365;
yrange=15:82;
dx=abs(Ds_low(1,2)-Ds_low(1,1)).*1e3;
dz=abs(De_low(2,1)-De_low(1,1));
%%%%%假设混合层内地转流均一
ualong2=repmat(ualong1,[size(Ds_low,1) 1]);
u=ualong2(yrange,xrange);
Ds_low=Ds_low(yrange,xrange);
De_low=De_low(yrange,xrange);
PotDen_low=PotDen_low(yrange,xrange);

rho1=PotDen_low+1000;
b=g.*(1-rho1./1025);
dbdz1=(b(2:end,:)-b(1:end-1,:))./dz;
N2=v2rho_2d(dbdz1);

bx=u2rho_2d((b(:,2:end)-b(:,1:end-1))./dx);
Ri=(N2.*f.^2)./(bx.^2);

Nthreshold=min(abs(N2(:)));
N2(N2==0)=Nthreshold;
N2(N2<0)=Nthreshold;
wqg=alongtrack_QGomega_MVP(u,b,N2,f,dx,dz,10000,1e-8);

% a=b(yrange,xrange);
a=u;
mean(a(:))
% contourf(Ds_low,De_low,Q)
% contourf(Ds_low,De_low,PotDen_low)

% x1=Ds_low;
% z1=De_low;

% w1=wqg;
% N1=N2;
% den1=PotDen_low;
% load mvppos.mat

target1=subplot(2,3,1)
pcolor(Ds_low,De_low,wqg);shading flat;colorbar;hold on
colortable = textread('MPL_seismic.txt');
% colortable = textread('MPL_viridis.txt');
contour(Ds_low,De_low,PotDen_low,[23.7,23.9,24],'color','k','linewi',2);
colormap(target1,colortable)
caxis([-2e-3 2e-3])
title('w');

target2=subplot(2,3,2)
pcolor(Ds_low,De_low,N2);shading flat;colorbar;hold on
% colortable = textread('GMT_polar.txt');
colortable = textread('MPL_viridis.txt');
contour(Ds_low,De_low,PotDen_low,[23.7,23.9,24],'color','w','linewi',2);
colormap(target2,colortable)
caxis([0 5e-4])
title('N2');

target3=subplot(2,3,3)

pcolor(Ds_low,De_low,PotDen_low);shading flat;colorbar;hold on
% colortable = textread('GMT_polar.txt');
colortable = textread('matlab_jet.txt');
contour(Ds_low,De_low,PotDen_low,[23.7,23.9,24],'color','w','linewi',2);
colormap(target3,colortable)
caxis([23.6 24.2])
title('rho');


%%
% bx=u2rho_2d((b(:,2:end)-b(:,1:end-1))./dx);
% Ri=(N2.*f.^2)./(bx.^2);
Ri(Ri<0)=0.25;
target4=subplot(2,3,4)

Ri1=Ri;
Ri1(Ri<=0)=0.25;Ri1(Ri<=0.25&&Ri>0)=0.5;Ri1(Ri>0.25&&Ri<=0.75)=0.75;Ri1(Ri>0.75)=1;
map = [1 0 0; 1 1 0; 0 0 1; 0 1 0];
% pcolor(Ds_low,De_low,log(Ri));shading interp;colorbar;hold on
pcolor(Ds_low,De_low,Ri1);shading interp;colorbar;hold on

% colortable = textread('GMT_polar.txt');
colortable = textread('matlab_hsv.txt');
contour(Ds_low,De_low,PotDen_low,[23.7,23.9,24],'color','k','linewi',2);
contour(Ds_low,De_low,log(Ri),[log(0.25) log(0.25)],'color','w','linewi',1);
colormap(target4,colortable)
colormap(target4,flipud(map))
% caxis([log(0.24) 2])
caxis([0 1])
title('Ri');

scale=Ri.^(-0.5);
target5=subplot(2,3,5)

pcolor(Ds_low,De_low,scale);shading interp;colorbar;hold on
% colortable = textread('GMT_polar.txt');
colortable = textread('MPL_ocean.txt');
contour(Ds_low,De_low,PotDen_low,[23.7,23.9,24],'color','k','linewi',2);
colormap(target5,flipud(colortable))
caxis([0 1])
title('wiso/wqg');

% target6=subplot(2,3,6)
% pcolor(Ds_low,De_low,wqg.*b);shading interp;colorbar;hold on
% % colortable = textread('GMT_polar.txt');
% colortable = textread('MPL_seismic.txt');
% contour(Ds_low,De_low,PotDen_low,[23.7,23.9,24],'color','k','linewi',2);
% colormap(target6,(colortable))
% caxis([-1e-5 1e-5])
% title('buoyancy flux');
% BF=wqg.*b;BF1=mean(BF,1);BF2=mean(BF1,2);

[N,L]=size(wqg);

for layer=1:N
     what=fftshift(fft((wqg(layer,:)),L));
     bhat=fftshift(fft(b(layer,:),L));
     bf(layer,:)=1./(length(b(layer,:)).^2).*real(what.*conj(bhat));
     kx=2.*pi./dx;
     kx1=[0:L/2-1 -L/2:-1]./L.*kx;
end
kx2=kx1(1,1:(L+1)/2);
bf1=bf(:,(L+1)/2:end);
z1=De_low(:,1);
z2=repmat(z1,[1,size(bf1,2)]);
kx3=repmat(kx2,[size(z2,1) 1]);

target6=subplot(2,3,6)
pcolor(kx3(:,2:end),z2(:,2:end),bf1(:,2:end));shading interp;colorbar;hold on
% colortable = textread('GMT_polar.txt');
colortable = textread('MPL_RdYlGn.txt');
colormap(target6,flipud(colortable))
caxis([-1e-11 1e-11])
title('buoyancy flux cross-spectrum');
saveas(gcf,'taiwanmvp2','png')
% bf1zavg=mean(bf1,1);
% plot(kx2,bf1zavg);
% set(gca,'Xscale','log')
