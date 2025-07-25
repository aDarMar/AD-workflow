close all; clear; clc;

global main_fold

main_fold = cd;
addpath('functions')
addpath('functions\Wing_Design');
addpath('functions\Wing_Design\Wing_Functions');
addpath('functions\Weight_estimation');
addpath('functions\Sizing_Plot'); 
addpath('functions\Classes'); addpath('functions\Wing_Design\Wing_Functions'); addpath('functions\Wing_Design\Wing_Functions\Grafici');
%% Lettura Nomi Aerei da File

name_list = dir('statistical_data/aircrafts');

nAero = length( name_list ) - 2;
Airl = AirData_class.empty;
for iAero = 1:nAero
    Airl(iAero) = AirData_class( [name_list(iAero+2).folder,'\',name_list(iAero+2).name] );
end
fun_plot_air( Airl,nAero );
[a,b]   = linear_regressions(Airl,nAero);
reg_cfs = [a,b];
fig_W = MTOMvsEM_plot(Airl,nAero,[a,b],[0.0810,1.0730;0.0913,1.0425]);
%% Aircraft Class & TLARs Definition
TLARs_path = [main_fold,'\tlars\TLARs.txt'];
Des_Air = Air_Design( TLARs_path,main_fold );

%% Weights Estimation
[Des_Air.MTOM_est0,Des_Air.EM_est0] = Class0_MTOM( Des_Air,reg_cfs,fig_W );

%% Sizing
% Input Data
CLmax_TO_vett = [ 2, 2.1, 2.2 ];    sigma = 1;
CLmax_CR_vett = [ 1.4,1.5,1.6 ];
CLmax_LND_vett = [ 2.1, 2.3, 2.5 ]; sigma = 1;

TisaoT50 = 1/0.8; phi_v = [1,1,1,0.85];
V_cr_vet = [Des_Air.TLARs.cruise.V, 1,Des_Air.TLARs.cruise.V,236];
h_cr_vet = [39000/3.3,Des_Air.TLARs.cruise.h, Des_Air.TLARs.cruise.h,11277];
n_fact   = [ 1,1.3, 1.3,1];
% Initialization
iS = 1;
Des_Air.SizHis(iS).WoS = 550;   % First Guess WoS [Kg/m^2]
Des_Air.SizHis(iS).S   = Des_Air.MTOM_est0/Des_Air.SizHis(iS).WoS;
%
[~, ~, ~, rho_h] = atmosisa( h_cr_vet(2) ); V_cr_vet(2) = sqrt( 2*Des_Air.SizHis(iS).WoS*9.81/(rho_h*2.1) )*1.2; % VFTO_calc
%
[Des_Air.SizHis(iS).CD0,Des_Air.SizHis(iS).Swet]   = Des_Air.polar_est( Des_Air.SizHis(iS).S,Des_Air.MTOM_est0 );
fig_ri       = figure('Name','Sizing Plot');

fig_aux = figure('Name','Auxiliary Figure for Sizing');
[ sizPLT_ax,ch_idxs,RoC_vt ] = sizing_plot(Des_Air,iS,CLmax_TO_vett,...
    CLmax_LND_vett,CLmax_CR_vett,sigma,TisaoT50,...
    V_cr_vet,h_cr_vet,n_fact,phi_v,fig_ri,fig_aux);

iS = 2;

flag = 1; tol = 1e-2;


while flag
    if iS > 2
        %
        [~, ~, ~, rho_h] = atmosisa( h_cr_vet(2) ); V_cr_vet(2) = sqrt( 2*Des_Air.SizHis(iS-1).WoS*9.81/(rho_h*2.1) )*1.2; % VFTO_calc [Optional]
        %
        [ sizPLT_ax,ch_idxs,RoC_vt ] = sizing_plot( Des_Air,iS-1,CLmax_TO_vett,...
            CLmax_LND_vett,CLmax_CR_vett,sigma,TisaoT50,...
            V_cr_vet,h_cr_vet,n_fact,phi_v,fig_ri,fig_aux,RoC_vt,ch_idxs );
    end
    % Plots a line corresponding to the assumed WoS
    lin           = plot( sizPLT_ax,Des_Air.SizHis(iS-1).WoS*[1,1],[0,1] );
    lin.LineStyle = '-'; lin.LineWidth = 1.0; lin.DisplayName = ['W/S = ',num2str(Des_Air.SizHis(iS-1).WoS)];
    col           = rand(1,3); lin.Color = col;
    legend( sizPLT_ax,'Interpreter','Latex' );
    
    tmp = input('Choose W/T');
    if tmp == -1 && iS > 2
        Des_Air.SizHis(iS-1).ToW = Des_Air.SizHis(iS-2).ToW;
    else
        Des_Air.SizHis(iS-1).ToW = tmp;
    end
    lin_pt = plot( Des_Air.SizHis(iS-1).WoS,Des_Air.SizHis(iS-1).ToW );
    lin_pt.LineStyle       = 'none';  lin_pt.Marker = 'o'; lin_pt.MarkerSize = 6;
    lin_pt.MarkerEdgeColor = [0, 1, 1];
    lin_pt.DisplayName     = 'Sizing Point';

    tmp = input('Choose W/S');
    if tmp == -1
        Des_Air.SizHis(iS).WoS = Des_Air.SizHis(iS-1).WoS;
    else
        Des_Air.SizHis(iS).WoS = tmp;
    end
    % Prelimianry Drag and S Update
    Des_Air.SizHis(iS).S                             = Des_Air.MTOM_est0/Des_Air.SizHis(iS).WoS;
    [Des_Air.SizHis(iS).CD0,Des_Air.SizHis(iS).Swet] = Des_Air.polar_est( Des_Air.SizHis(iS).S,Des_Air.MTOM_est0 ); %polar_est(Sizing(iS).S,MTOM_it0);
    
    err.WoS = abs( (Des_Air.SizHis(iS).WoS - Des_Air.SizHis(iS-1).WoS)/Des_Air.SizHis(iS-1).WoS ) ;
    %abs.WoT = ( Sizing(iS).WoT - Sizing(iS-1).WoT ) ;
    flag = err.WoS > tol; %&& abs.WoT <tol;
    iS = iS + 1;
    %hold off
end
sizPLT_ax_cf = sizing_plot_cfr( Airl,nAero,sizPLT_ax );
xlabel(sizPLT_ax,'W/S [Kg/m$^2$]','Interpreter','Latex','FontSize',16);
ylabel(sizPLT_ax,'T/W','Interpreter','Latex','FontSize',16)
%% Wing Design
Des_Air = Des_Air.final_out( ch_idxs,CLmax_TO_vett,CLmax_CR_vett,...
    CLmax_LND_vett,V_cr_vet,h_cr_vet );
fig_p = Des_Air.preliminary_polar_plot;
Wing_Design( Des_Air );
