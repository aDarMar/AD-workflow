function [MAC,k1,k2,k3,yroot,ykink,ytip,croot,ckink,ctip,xLE_root,xLE_kink, ...
    xLE_tip] = Wing_Planform(Sw,ARw,Mcruise,sweep,tr1,...
    ykink_b,tr2,eps_root,eps_tip,A1_Sw,Wi_cruise, ...
    Wf_cruise,hcruise,CL_cruise,CL_max,SM,ka,tc_root,tc_kink,tc_tip)
%INPUT
% Sw        : superficie alare
% ARw       : aspect ratio dell'ala
% Mcruise   : Mach di crociera
% sweep     : angolo di freccia
% tr1       : taper ratio c_kink/c_root (panel 1)
% ykink_b   : posizione del kink lungo lo span
% tr2       : taper ratio c_tip/c_root (panel 2)
% eps_root  : twist del profilo di root
% eps_tip   : twist del profilo di tip
% A1_Sw     : percentuatel dl'area del panel 1 rispetto alla totale
% Wi_cruise : peso inizio crociera
% Wf_cruise : peso fine crociera
% hcruise   : quota di crociera
% CL_cruise : CL in cruise
% CL_max    : CL massimo dell'ala
% SM        : coefficiente di sicurezza (0.05)
% ka        : coefficiente per il calcolo del t/c (dipende dal profilo)
% tc_root   : thickness ratio del profilo di root
% tc_kink   : thickness ratio del profilo di kink
% tc_tip    : thickness ratio del profilo di tip
%OUTPUT
% MAC       : corda media aerodinamica
% k1        : coeffiente per il triangolo 1 dell'ala
% k2        : coeffiente per il triangolo 2 dell'ala
% k3        : coeffiente per il triangolo 3 dell'ala
% yroot     : posizione della root lungo lo span (0)
% ykink     : posizione della kink lungo lo span
% ytip      : posizione della tip  lungo lo span
% croot     : corda del profilo di root
% ckink     : corda del profilo di kink
% ctip      : corda del profilo di tip
% xLE_root  : posizione del LE della root
% xLE_kink  : posizione del LE della kink
% xLE_tip   : posizione del LE della tip

bw       = sqrt(Sw*ARw);
yroot    = 0;    
ykink    = ykink_b*bw/2;      
ytip     = bw/2;
bkink    = bw*ykink_b/2;        
A1d      = A1_Sw*Sw;
A        = [bkink bkink;tr1 -1];  
noti     = [2*A1d ; 0];     
Ainv     = inv(A);
croot    = Ainv(1,:)*noti;     
ckink    = Ainv(2,:)*noti;     
ctip     = tr2*croot;
xLE_root = 0;                                xTE_root = xLE_root + croot;
xLE_kink = xLE_root+tand(sweep)*ykink;       xTE_kink = xLE_kink + ckink;
xLE_tip  = xLE_root+tand(sweep)*ytip;        xTE_tip  = xLE_tip + ctip;

% Wing Aerea Check
A1   = .5*(ckink+croot)*(ykink-yroot);
A2   = .5*(ckink+ctip)*(ytip-ykink);
Atot = 2*(A1+A2);   diff = Sw-Atot;
fprintf("La differenza tra Sw scelta nel Sizing e la Sw calcolata con il" + ...
    "metodo delle aeree è: %.3f\n",diff);

% Equivalent Wing
croot_eq = ((((Sw/(0.5*bw*(1-yroot/bw)))-2*ctip)/(0.5*bw*(1-yroot/bw)))*0.5* ...
    bw*yroot/bw)+(Sw/(0.5*bw*(1-yroot/bw)))-ctip;
ctip_eq     = ctip;
xTE_root_eq = xLE_root + croot_eq;
xTE_tip_eq  = xLE_tip + ctip_eq;
tr_eq       = ctip_eq/croot_eq;

sweep_eq    = atan((xLE_tip-xLE_root)/(ytip-yroot))*57.3;
sweep_eq_c4 = atan(tand(sweep_eq)-(4/ARw)*(.25*(1-tr_eq)/(1+tr_eq)))*57.3;
sweep_eq_c2 = atan(tand(sweep_eq)-(4/ARw)*(.5*(1-tr_eq)/(1+tr_eq)))*57.3;

slop_c      = (ctip_eq-croot_eq)*2/bw;
slop_tw     = (eps_tip-eps_root)/ytip;
slop_xle    = (xLE_tip-xLE_root)/ytip;

% Weighted Average Area
S1 = .5*croot*(ykink-yroot);
S2 = .5*ckink*(ykink-yroot)+.5*ckink*(ytip-ykink);
S3 = .5*ctip*(ytip-ykink);
Stot = 2*(S1+S2+S3);
k1 = 2*S1/Stot;  k2 = 2*S2/Stot;  k3 = 2*S3/Stot;  ktot = k1 + k2 + k3;

% Panels and MAC
xLE_reg_1 = polyfit([yroot;ykink],[xLE_root;xLE_kink],1);
xLE_reg_2 = polyfit([ykink;ytip],[xLE_kink,xLE_tip],1);
AR_1      = ((ykink-yroot)^2)/A1;   
AR_2      = ((ytip-ykink)^2)/A2; 
sweep_c4_1 = 57.3*atan(tand(sweep_eq)-(4/AR_1)*(.25*(1-(ckink/croot))/ ...
    (1+(ckink/croot))));
sweep_c4_2 = 57.3*atan(tand(sweep_eq)-(4/AR_2)*(.25*(1-(ctip/ckink))/ ...
    (1+(ctip/ckink))));
MAC  = (2/3)*(croot_eq+ctip_eq-(croot_eq*ctip_eq/(croot_eq+ctip_eq)));
YMAC = (1/3)*((1+2*tr_eq)/(1+tr_eq))*ytip;
XMAC = YMAC*tand(sweep_eq);
if YMAC <= ykink
    xLE_MAC = xLE_reg_1(2)+xLE_reg_1(1)*YMAC;
else
    xLE_MAC = xLE_reg_2(2)+xLE_reg_2(1)*YMAC;
end

%% Airfoil Selection
Wav_cruise   = (Wi_cruise+Wf_cruise)/2;
CL_cruise_av = CL_cruise/.95;
CL_max_W     = CL_max/.95;
[~, a]       = atmosisa(hcruise);
vcruise      = a*Mcruise;

cl_meanairf_cruise = CL_cruise_av/.9;   cl_max_meanairf = CL_max_W/.9;
M_DD          = Mcruise+SM;
tc_mean       = (-cl_meanairf_cruise-10*M_DD*cosd(sweep_eq_c4)^3+10*ka* ...
    cosd(sweep_eq_c4)^2)/(10*cosd(sweep_eq_c4));
M_design      = Mcruise*cosd(sweep);
tc_mean_input = tc_root*k1+tc_kink*k2+tc_tip*k3;
if tc_mean_input>tc_mean
    fprintf("Errore! La distribuzione degli spessori non rispetta la condizione..." + ...
        "limite. Correggi")
end
end