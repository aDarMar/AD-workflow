function [dCd0] = dCd0_flaps( flap_obj,deltaF,wing_obj )
%DCD0_FLAPS Summary of this function goes here
%   INPUT: 
%       cfoc:
%       

load('/Grafici/dCd0_HL_slotted.mat');
%% Delta 1
del1    = scatteredInterpolant( tc_d1,cf_c_d1,d1 );
tc_avg  = 0.5*( flap_obj.root.tc + flap_obj.tip.tc ); cfoc_avg = 0.5*( flap_obj.root.cfoc + flap_obj.tip.cfoc );
delta_1 = del1( tc_avg,cfoc_avg );
%% Delta 2
del2    = scatteredInterpolant( t_c_delta2,df_d2,d2 );
delta_2 = del2( tc_avg,deltaF );
%% Delta 3
del3    = scatteredInterpolant( taper_d3,bf_b_d3,d3 );
 [croot_eq,ctip_eq,~] = Equivalent_Wing( wing_obj.Sw,wing_obj.ARw,...
     wing_obj.panels(1).root.eps,wing_obj.panels(end).tip.eps );
 TR_e = ctip_eq/croot_eq;
delta_3 = del3( TR_e,2*(flap_obj.tip.yglob-flap_obj.root.yglob)/wing_obj.bw );
%% dCd0
dCd0 = delta_1*delta_2*delta_3;
end

