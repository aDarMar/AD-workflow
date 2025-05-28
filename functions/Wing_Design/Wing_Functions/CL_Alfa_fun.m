function [CL_Alfa] = CL_Alfa_fun(Mean_Cl_alfa,Freccia_LE,AR,Freccia_c2,Mach)
% froma2AFun: funzione che corregge il Cl_alpha del profilo medio nel
% CL_alpha dell'ala
% Input: a, sweep_le, AR, sweep_c2, M. Gli angoli sono in gradi

if Freccia_LE>10
    k = Mean_Cl_alfa*57.3*cos(Freccia_c2/57.3); % ao*cos(Sweep) <- converte il risultato in rad^-1
    j = 1 - (Mach^2)*(cos(Freccia_c2/57.2))^2;  % 1-M^2cos^2(sweep)

    CL_Alfa = k/( sqrt(j+ ( k/(pi*AR) )^2) +  k/(pi*AR) );  CL_Alfa  =  CL_Alfa /57.3;

else
    
    CL_Alfa=(Mean_Cl_alfa*57.3)/(1+((57.3*Mean_Cl_alfa)/(pi*AR))/57.3);
    
end
end
