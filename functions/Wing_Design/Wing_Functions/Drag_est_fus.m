function [CD_tot_f,CD_0_n] = Drag_est_fus(M,h,Lf,df,Type_of_surface_f,db,...
    Dcd_windshield,h_l_075,S_wet_F,Sw,Ln,FF_n,SwetN)
%DragEstFusNac funzione che calcola il CD0 di fusoliera e delel nacelles.
%Nacelles
%  Ln: fan cowling length
%  FF_n: form factor nacelles
%  SwetN: area bagnata nacelles



[~,a,~,rho,ni,~] = atmosisa(h);

K_rough_f      = K_rough_calc(Type_of_surface_f);
k_rough_f      = K_rough_f*10^-3; % m

% Nacelles
Re_n           = Re_calc(M,a,Ln,ni);
Re_cut_off_n   = Re_cut_off_calc(Ln,k_rough_f,M);
cf_turb_n      = cf_turb_calc(Re_n,Re_cut_off_n,M);
CD_0_n         = cf_turb_n*FF_n*SwetN/Sw;

% Fusoliera
S_front_f = pi*(df*0.5)^2;
FF_f           = 1 + ( (60)/(Lf/df)^3 ) + ( Lf/df )/400;
Re_f           = Re_calc(M,a,Lf,ni);
Re_cut_off_f   = Re_cut_off_calc(Lf,k_rough_f,M);
cf_turb_f      = cf_turb_calc(Re_f,Re_cut_off_f,M);
CD_0_f         = cf_turb_f*FF_f*S_wet_F/Sw;
CD_upsweep     = 0.075*(h_l_075)*pi*( (df/2)^2 )/Sw;
CD_base        = ( S_front_f/Sw )*( (0.029*(db/df)^3 )/( CD_0_f*Sw/S_front_f )^0.5 );
CD_windshield  = Dcd_windshield*S_front_f/Sw;
CD_tot_f = CD_0_f+CD_upsweep+CD_base+CD_windshield;
end