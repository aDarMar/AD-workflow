function [toc] = Airfoil_Selection( aero_des,SM,ka,toc,c,y)
%AIRFOIL_SELECTION Function that gives the thickness distribution for a
%given Mach number
%INPUT
%   aero_des: Air_Design class object containing the equivalent wing class
%       object

Mcruise   = aero_des.TLARs.cruise.M;
CL_cruise = aero_des.CL_cr; % CL of the A/C
CL_max    = aero_des.CLmax_cr;
sweep     = aero_des.sweepw;
yroot = 0;    ykink = y(2); ytip = y(3);
croot = c(1); ckink = c(2); ctip = c(3);

sweep_c4 = aero_des.equiv_wing.sweepChange( ...
    aero_des.equiv_wing.panels.sweep,0,0.25,aero_des.equiv_wing.AR,...
    aero_des.equiv_wing.TR );
%% Weighted Average Area
S1 = .5*croot*(ykink-yroot);
S2 = .5*ckink*(ykink-yroot)+.5*ckink*(ytip-ykink);
S3 = .5*ctip*(ytip-ykink);
Stot = 2*(S1+S2+S3);
k1 = 2*S1/Stot;  k2 = 2*S2/Stot;  k3 = 2*S3/Stot;  ktot = k1 + k2 + k3;

%% CL cruise
CL_cruise_av = CL_cruise/.95; % CL cruise of the isolated Wing
CL_max_W     = CL_max/.95;

%% t/c mean
% Cl 2D from wing CL 3D
cl_meanairf_cruise = CL_cruise_av/.9;   cl_max_meanairf = CL_max_W/.9;
M_DD          = Mcruise+SM;
tc_mean       = ( -cl_meanairf_cruise-10*M_DD*cosd(sweep_c4)^3+10*ka* ...
    cosd(sweep_c4)^2)/( 10*cosd(sweep_c4) );
M_design      = Mcruise*cosd(sweep);
err_tc = 1;%while flag
% Input Cycle for t/c
while err_tc
    while 1 > 0
        disp('------------------------------')
        disp('Root Kink Tip');
        disp(toc(:)');
        disp('Choose section to modify. Press any other number to skip')
        sec_CHS = input('>>');
        switch sec_CHS
            case 1
                disp('Modifyng Root');
                disp('Insert t/c')
                temp   = input('>>');
                toc(1) = temp;
            case 2
                disp('Modifyng Kink');
                disp('Insert t/c')
                temp   = input('>>');
                toc(2) = temp;
            case 3
                disp('Modifyng tip');
                disp('Insert t/c')
                temp   = input('>>');
                toc(3) = temp;
            otherwise
                break
        end
    end
    tc_mean_input = toc(1)*k1+toc(2)*k2+toc(3)*k3;
    err_tc = 0;
    disp(['t/c required for compressibility: ',num2str(tc_mean)])
    disp(['t/c assigned: ',num2str(tc_mean_input)])
    if tc_mean_input>tc_mean
        fprintf("Errore! La distribuzione degli spessori non rispetta la condizione..." + ...
            "limite. Correggi")
        err_tc = 1;
    else
    
    disp('Are you satisfied?')
    disp('1.Yes')
    chs = input('>>');
    if chs == 1
        err_tc = 0;
    else
        err_tc = 1;
    end
    end
end
end

