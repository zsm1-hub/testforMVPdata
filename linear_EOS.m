function linear_EOS(temp1,salt1)
rho0=1025;T0=10;S0=35;TCOFF=1.7e-4;SCOEF=7.6e-4;
rho1=rho0.*(1-TCOFF.*(temp1-T0)+SCOEF.*(salt1-S0));
return