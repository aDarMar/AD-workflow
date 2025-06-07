function CDlow = CDlow_Mach(a,cd_mean,a_w,CLw_s,TR,AR,eps_ae,Cl_a)
% CDSTALL Calcola il drag totale a basso Mach (senza drag d’onda)
% 
% INPUT:
%   a        - Vettore angoli d'attacco per il profilo 2D
%   cd_mean  - Vettore corrispondente dei cd 2D
%   a_w      - Alpha dell’ala da considerare
%   CLw_s    - Portanza dell’ala in low Mach
%   AR       - Aspect Ratio
%   TR       - Taper Ratio
%   eps_ae   - Deviazione incidenza effettiva
%   Cl_a     - Derivata Cl media (2D)
%
% OUTPUT:
%   CDstall  - Coefficiente di drag totale in condizioni di stallo
erroreMAX = 0.05; % imposto il minimo errore
for n=1:5
    p     = polyfit(a,cd_mean,n);
    cdfit = polyval(p,a_w); % cd del profilo 2d ottenuto con il polinomio 
    err   = norm(cd_mean-cdfit);
    if err<erroreMAX
        break;
    end
end
u = interpolateFromCSV('AR*.csv', TR, AR);
v = interpolateFromCSV('TR*.csv', AR, TR);
w = interpolateFromCSV('ctcr*.csv', AR, TR);


t1 = CLw_s^2/(pi*AR*u);%c'è un fattore s che non sappiamo cosa significa anche nell'excel non viene calcolato
t2 = v*CLw_s*eps_ae*Cl_a;
t3 = (eps_ae+Cl_a)^2*w;

CDi = t1+t2+t3; 

CDlow = cdfit+ CDi; %non c'è contributo di wave perchè siamo a basso mach
end