function [lon,lat,ugos,vgos,sla]=read_AVISO(fname)
% fname='cmems_mod_glo_phy_my_0.083deg_P1D-m_1716049535443.nc';
nc=netcdf(fname,'r');
ncdisp(fname,'/','full')
lon=nc{'longitude'}(:);
lat=nc{'latitude'}(:);
ugos=nc{'ugos'}(:);
vgos=nc{'vgos'}(:);
sla=nc{'sla'}(:);