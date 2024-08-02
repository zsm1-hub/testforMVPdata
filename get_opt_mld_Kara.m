function [ILD,MLD]=get_opt_mld_Kara();
%%%A. Birol Kara, 2000, JGR, optiminal definition for MLD

deltaT=0.2;deltaT2=deltaT.*0.1;

for ii=1:size(z2,2)
    T2=temp2(:,ii);
    zzz=z2(:,ii);
    mask1=~isnan(T2);
    T3=T2(mask1);
    Z3=zzz(mask1);

    a1=abs(abs(Z3)-10);


    zpos1=find(a1==min(a1));
    %%%step 1, choose 10m as a reference level
    Tref=T3(zpos1);

    %%% step 2, search uniform temperature region
    a2=abs(T3-Tref)<deltaT2;
    a3=find(a2==1);
    a3(a3<=10)=[];
    zpos2=min(a3);
    Tref=T3(zpos2);

    if isempty(zpos2)==1

        ILD(:,ii)=NaN;
    else
        % a22=~isnan(T2);
        %% step 3, update point?step 4,critera
        if T2(zpos2)>T2(zpos2+1)
            Tb=Tref+deltaT;
            ILD(:,ii)=interp1(T3(zpos2:end),Z3(zpos2:end),Tb);
        else
            Tb=Tref-deltaT;
            ILD(:,ii)=interp1(T3(1:zpos2),Z3(1:zpos2),Tb);
        end
        % disp()
        disp(ILD)
        % ILD(:,ii)=interp1(T2(a22),zzz(a22),Tb);

    end
end


return