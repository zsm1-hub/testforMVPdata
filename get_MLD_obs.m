function [MLDmix,MLDt,MLDr]=get_MLD_obs(temp,rho,z)
for ii=1:size(temp,2)
    a1=temp(:,ii);
    b1=rho(:,ii);
    zdot=min(find(~isnan(a1)==1));
    sur_temp=a1(zdot);
    sur_rho=b1(zdot);
    MLD_temp=sur_temp-0.2;
    MLD_rho=sur_rho+0.03;
    zzz=z(:,ii);
    zdot1=~isnan(a1);
    a2=a1(zdot1);
    b2=b1(zdot1);
    zzz2=zzz(zdot1);
    if isempty(a2)==1
        MLDt(:,ii)=nan;
        MLDr(:,ii)=nan;
    else
        MLDt(:,ii)=interp1(a2,zzz2,MLD_temp,'nearest');
        MLDr(:,ii)=interp1(b2,zzz2,MLD_rho,'nearest');
    end

end
    for ii=1:size(temp,2)
    MLDmix(:,ii)=min([MLDt(:,ii),MLDr(:,ii)]);
    end
return