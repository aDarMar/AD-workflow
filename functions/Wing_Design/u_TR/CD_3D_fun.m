function CDfun = CDw(CLw_c, Cd, ka, sweep_c4, tc_mean,M_cruise,TR,AR, eps_ae,Cl_alpha_m)
% CDw Calcola il coefficiente di drag totale di un'ala in cruise
% 
% INPUT:
%   CLw_c       - Coefficiente di portanza in crociera
%   Cd          - Drag del profilo medio (2D)
%   ka          - Correttivo per profili supercritici
%   sweep_c4    - Freccia alla c/4 [rad]
%   tc_mean     - Spessore/corda medio
%   M_cruise    - Mach in crociera
%   TR          - Taper ratio
%   AR          - Aspect ratio
%   eps_ae      - Deviazione incidenza effettiva
%   Cl_alpha_m  - Derivata portanza media
%
% OUTPUT:
%   CDfun       - Drag totale


MDD = ka/(cos(sweep_c4)) - (tc_mean/(cos(sweep_c4)^2))-(CLw_c/(10*(cos(sweep_c4))^3));
Mcrit = MDD-(0.1/80)^(1/3);

DeltaCd_wave = 20*(M_cruise-Mcrit)^4;

%calcolo dei coefficienti u,v,z per il CDi con la function
%interpolatefromcsv fatta a parte
u = interpolateFromCSV('AR*.csv', TR, AR);
v = interpolateFromCSV('TR*.csv', AR, TR);
w = interpolateFromCSV('ctcr*.csv', AR, TR);

t1 = CLw_c.^2/(pi*AR*u);%c'è un fattore s che non sappiamo cosa significa anche nell'excel non viene calcolato
t2 = v*CLw_c*eps_ae*Cl_alpha_m;
t3 = (eps_ae+Cl_alpha_m)^2*w;

CDi = t1+t2+t3;

CDfun = Cd + CDi + DeltaCd_wave;
end