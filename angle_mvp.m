function [ualong,uacross,angle]=angle_mvp(lon,lat,u,v)
distx=spheric_dist(lat(2:end),lat(2:end),lon(2:end),lon(1:end-1));
disty=spheric_dist(lat(2:end),lat(1:end-1),lon(2:end),lon(2:end));
angle=atand(disty./distx);
angle(isnan(angle))=nanmean(angle);
angle(angle>90)=nanmean(angle);

u(isnan(u))=nanmean(u);
v(isnan(v))=nanmean(v);

u1=0.5*(u(2:end)+u(1:end-1));
v1=0.5*(v(2:end)+v(1:end-1));

ux=cosd(angle).*u1-sind(angle).*v1;
uy=sind(angle).*u1+cosd(angle).*v1;

ualong=v2rho_2d(ux);
uacross=v2rho_2d(uy);
return

