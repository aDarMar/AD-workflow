function [croot_eq,ctip_eq,sweep_eq] = Equivalent_Wing( Sw,ARw,eps_root,eps_tip )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

%% Equivalent Wing
croot_eq = ((((Sw/(0.5*bw*(1-yroot/bw)))-2*ctip)/(0.5*bw*(1-yroot/bw)))*0.5* ...
    bw*yroot/bw)+(Sw/(0.5*bw*(1-yroot/bw)))-ctip;
ctip_eq     = ctip;
xTE_root_eq = xLE_root + croot_eq;
xTE_tip_eq  = xLE_tip + ctip_eq;
tr_eq       = ctip_eq/croot_eq;

sweep_eq    = atan((xLE_tip-xLE_root)/(ytip-yroot))*57.3;
% sweep_eq_c4 = atan(tand(sweep_eq)-(4/ARw)*(.25*(1-tr_eq)/(1+tr_eq)))*57.3;
% sweep_eq_c2 = atan(tand(sweep_eq)-(4/ARw)*(.5*(1-tr_eq)/(1+tr_eq)))*57.3;

slop_c      = (ctip_eq-croot_eq)*2/bw;
slop_tw     = (eps_tip-eps_root)/ytip;
slop_xle    = (xLE_tip-xLE_root)/ytip;

end