function [mld]=get_mld_taiwanPIO(rho2,z2);
%%%A. Birol Kara, 2000, JGR, optiminal definition for MLD

% deltaT=0.2;deltaT2=deltaT.*0.1;
%%%for taiwan PIO
deltaR=0.125

for ii=1:size(z2,2)
    R2=rho2(:,ii);
    zzz=z2(:,ii);
    mask1=~isnan(R2);
    R3=R2(mask1);
    Z3=zzz(mask1);

    a1=abs(abs(Z3)-5);


    zpos1=find(a1==min(a1));
    %%%step 1, choose 10m as a reference level
    Rref=max(R3(zpos1));
    % Rref=min(R3(zpos1));

    %%% step 2, exceed 0.125
    Rmix=Rref+deltaR;
    if isempty(Z3)==1
        mld(ii)=NaN;
    else
        if Rmix>max(R3)
            mld(ii)=min(Z3);
        else
            mld(ii)=min(interp1(R3,Z3,Rmix));
        end
    end
end


return