function [GIm,GISIm,SIm,ISIm,Stablem]=get_thomas_meassure(Rib,xig,f);
GIm=Rib;GIm(Rib<-1)=1;GIm(GIm~=1)=nan;
GISIm=(double(-1<Rib).*double(Rib<0));GISIm(GISIm~=1)=nan;
SIm=double(xig<0).*double(Rib>0).*double(Rib<1)+...
    double(xig>0).*double(Rib>0).*double(Rib<(f./xig));SIm(SIm~=1)=nan;
ISIm=double(Rib>1).*double(Rib<(f./xig));ISIm(ISIm~=1)=nan;
Stablem=double(Rib>(f./xig));

% subplot(211)
% hold on;
% scatter(x2.*GIm,z2.*GIm,[],[0, 0.5, 0],"filled")
% scatter(x2.*GISIm,z2.*GISIm,[],[0, 0, 0.5],"filled")
% scatter(x2.*SIm,z2.*SIm,[],[0.8, 0.6, 0.4],"filled")
% scatter(x2.*ISIm,z2.*ISIm,[],[0.7, 0.7, 0.7],"filled")
% scatter(x2.*Stablem,z2.*Stablem,[],[1, 0.1, 0.1],"filled")
% 
% contour(x2,z2,filter_rho2,[22.5:pycnal:24],'linewi',.5,'linestyle','-','color',colorcon);
% contour(x2,z2,mask1,[1 1],'linewi',1.5,'linestyle','--','color',colorcon);
% plot(mld_CDx(1,:),mld,'color','r','linestyle','--','LineWidth',1.5)
% 
% 
% subplot(212)
% hold on;
% scatter(x2.*GIm,z2.*GIm,[],[0, 0.5, 0],"filled")
% scatter(x2.*GISIm,z2.*GISIm,[],[0, 0, 0.5],"filled")
% scatter(x2.*SIm,z2.*SIm,[],[0.8, 0.6, 0.4],"filled")
% scatter(x2.*ISIm,z2.*ISIm,[],[0.7, 0.7, 0.7],"filled")
% % scatter(x2.*Stablem,z2.*Stablem,[],[1, 0.1, 0.1],"filled")
% 
% contour(x2,z2,filter_rho2,[22.5:pycnal:24],'linewi',.5,'linestyle','-','color',colorcon);
% contour(x2,z2,mask1,[1 1],'linewi',1.5,'linestyle','--','color',colorcon);
% plot(mld_CDx(1,:),mld,'color','r','linestyle','--','LineWidth',1.5)
% 
return