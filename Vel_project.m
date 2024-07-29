function [ualong,uacross,angle]=Vel_project(lon,lat,u,v)
distx=spheric_dist(lat(2:end),lat(2:end),lon(2:end),lon(1:end-1));
disty=spheric_dist(lat(2:end),lat(1:end-1),lon(2:end),lon(2:end));

angle=atand(disty./distx);
u1=0.5*(u(2:end,:)+u(1:end-1,:));
v1=0.5*(v(2:end,:)+v(1:end-1,:));
ANGLE=repmat(angle,[1,size(u1,2)]);
ux=cosd(ANGLE).*u1-sind(ANGLE).*v1;
uy=sind(ANGLE).*u1+cosd(ANGLE).*v1;
ualong=v2rho_2d(ux);
uacross=v2rho_2d(uy);

return
