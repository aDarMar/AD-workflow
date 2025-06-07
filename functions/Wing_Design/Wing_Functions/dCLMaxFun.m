function [dCLmax] = dCLMaxFun(sweep, dy,M_inf,c)

load('Grafici\Fig5.18.mat');

if nargin ==4
    dy_c=dy/c;
else
    dy_c = dy;
end

Interpol= scatteredInterpolant(dY_Cscat,Mscat,Sweepscat,dCLmaxscat,'linear', 'boundary');
dCLmax=Interpol(dy_c, M_inf, sweep); 

