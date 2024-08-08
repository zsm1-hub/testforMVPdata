function [u,v]=wind_vec(ws,wd)
u=ws.*sin((wd+180)/360*2*pi);
v=ws.*cos((wd+180)/360*2*pi);
% 
% % tmp=270.0-atan2(v,u)*180.0/pi;
% % dir=mod(tmp,360.0);
% % ws=sqrt(u.^2+v.^2);