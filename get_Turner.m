function [R,Tu]=get_Turner(alpha,beta,T,S,X);
% R = alphaDT/betaDS
% Tu=arctanR
X=X.*1e3;
dT=(T(:,2:end)-T(:,1:end-1))./abs(X(:,2:end)-X(:,1:end-1));
dS=(S(:,2:end)-S(:,1:end-1))./abs(X(:,2:end)-X(:,1:end-1));
R=(alpha.*dT)./(beta.*dS);
Tu=atan(R);
return
