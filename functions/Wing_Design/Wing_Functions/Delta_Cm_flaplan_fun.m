function [Delta_Cm_flap_lan] = Delta_Cm_flaplan_fun(Delta_CL_max,cc,cfc,...
    Cl,S_wing,S_flap,AR,Delta_Cl_max,Freccia_c4, deltaF,bfob,TR)
%UNTITLED2 Summary of this function goes here
%   Delta_CL_max
%   cc = c'/c
%   cfc = cf/c medio
%   Cl: il CL max.
%   S_wing
%   S_flap
%   AR
%   Delta_Cl_max
%   Freccia_c4
%   deltaF: deflessione degli slat
if deltaF == 0
    Delta_Cm_flap_lan = 0;
else
    load("Grafici\mu1_tab.mat","cf_c","delta_f","mu_1");
    mu1 = interp2(cf_c,delta_f,mu_1,cfc/cc,deltaF);
    %u1 = interp2(cf_c,delta_f,mu_1,0.24898,40);
    load("Grafici\mu2_tab.mat","i");
    i.ExtrapolationMethod = 'boundary';
    mu2 = i(bfob,TR);
    load("Grafici\mu3_tab.mat","int");
    int.ExtrapolationMethod = 'boundary';
    mu3 = int(bfob,TR);


    Delta_Cm_flap_lan = mu2*(-mu1*Delta_CL_max*cc-(Cl+Delta_Cl_max*(1-S_flap/S_wing))*...
        (1/8)*cc*(cc-1))+0.7*AR/(1+AR/2)*mu3*Delta_Cl_max*tand(Freccia_c4);
end
end