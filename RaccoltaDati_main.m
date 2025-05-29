close all; clear; clc;

global main_fold

main_fold = cd;
addpath('functions')
addpath('functions\Wing_Design');
addpath('functions\Wing_Design\Wing_Functions');
addpath('functions\Weight_estimation');
addpath('functions\Sizing_Plot');
addpath('functions\Classes');
%% Lettura Nomi Aerei da File

name_list = dir('statistical_data/aircrafts');

nAero = length( name_list ) - 2;
Airl = AirData_class.empty;
for iAero = 1:nAero
    Airl(iAero) = AirData_class( [name_list(iAero+2).folder,'\',name_list(iAero+2).name] );
end

[a,b] = linear_regressions(Airl,nAero);
reg_cfs = [a,b];
MTOMvsEM_plot(Airl,nAero,[a,b],[0.0810,1.0730;0.0913,1.0425]);
%% Aircraft Class 
TLARs_path = [main_fold,'\tlars\TLARs.txt'];
Des_Air = Air_Design( TLARs_path,main_fold );
%% TLARs Definition
% % % Definizione TLARS ( da portare in una funzione a parte )
% % TLARS = read_TLARs( TLARs_path );

%% Weights Estimation
% Fuel fraction
% % [Mff,Mff_b,WLNDoWTO,WcroWTO] = fuel_fraction(TLARS);
% % 
% % M_pay  = TLARS.npax*215/2.2046; %[Kg]
% % M_crew = (TLARS.ncrew+TLARS.npil)*205/2.2046; %[Kg]
% % Mres = 0; Mfo = 0;
% % c = 1 - (1+Mres)*(1-Mff) - Mfo; d = M_pay + M_crew;

[Des_Air.MTOM_est0,Des_Air.EM_est0] = Class0_MTOM( Des_Air,reg_cfs );

% % Wmtom_reg   = (2:20)*1e4; %Wmtom_reg = Wmtom_reg*2.2046;                          % da [Kg] -> [lb]
% % Wempty_reg  = @(Wmtom_reg) ( 10.^( ( log10(Wmtom_reg*2.2046) - a )./b )/2.2046 ) ;      % in [kg]
% % Wempty_stat = @(Wmtom_reg) c*Wmtom_reg-d;
% % find_W      = @(Wmtom_reg) Wempty_reg(Wmtom_reg) - Wempty_stat(Wmtom_reg);
% % x0          = 120000; % [Kg]
% % MTOM_it0    = fzero( find_W,x0 );
% % % = Wempty_reg/2.2046; Wmtom_reg = Wmtom_reg/2.2046;       % da [lb] -> [Kg]
% % 
% % % DEBUGGGG 
% % %MTOM_it0 = 95000;
% % 
% % figure()
% % plot( Wmtom_reg,Wempty_reg( Wmtom_reg ),'--r' ); hold on
% % plot( Wmtom_reg,Wempty_stat( Wmtom_reg ),'k' );
% % 
% % % plot( Wmtom_reg,Wempty_reg,'--r' ); hold on
% % % plot( Wmtom_reg,c*Wmtom_reg-d,'k' );
% % 
% % axis equal
% % Polar Estimation
% WoS_it0 = 550; %[kg/m^2]
% S_it0 = MTOM_it0/WoS_it0;
% [CD0,Swet] = polar_est(S_it0,MTOM_it0);

%% Sizing
% Input Data
CLmax_TO_vett = [ 2, 2.1, 2.2 ];    sigma = 1;
CLmax_CR_vett = [ 1.4,1.5,1.6 ];
CLmax_LND_vett = [ 2.1, 2.3, 2.5 ]; sigma = 1;

TisaoT50 = 1/0.8; phi_v = [1,0.85];
V_cr_vet = [Des_Air.TLARs.cruise.V,236];
h_cr_vet = [Des_Air.TLARs.cruise.h,11277];

% Initialization
iS = 1;
Des_Air.SizHis(iS).WoS = 550;   % First Guess WoS [Kg/m^2]

Des_Air.SizHis(iS).S   = Des_Air.MTOM_est0/Des_Air.SizHis(iS).WoS;
% FINIREEEEE
[Des_Air.SizHis(iS).CD0,Des_Air.SizHis(iS).Swet]   = Des_Air.polar_est( Des_Air.SizHis(iS).S,Des_Air.MTOM_est0 );
fig_ri       = figure(); %axis([0,1000,0,1]); 
%ax_fig = axes('Parent', fig_ri); %axis([0,1000,0,1]); hold on;  
% grid minor; xlabel('W$_{TO}$ / S [Kg/m$^2$]','Interpreter','latex','FontSize',16);
% ylabel('T$_{TO}$ / W$_{TO}$ [-]','Interpreter','latex','FontSize',16);
%hold on;

fig_aux = figure();
[ch_idxs,ax_siz,Leg_siz,idx_chs] = sizing_plot(Des_Air.TLARs,CLmax_TO_vett,...
    CLmax_LND_vett,CLmax_CR_vett,Des_Air.MLndoMTo,sigma,...
    Des_Air.SizHis(iS).CD0,Des_Air.dCD0_wave,TisaoT50,Des_Air.MCroMTo,...
    V_cr_vet,h_cr_vet,phi_v,fig_ri,fig_aux);

iS = 2;

flag = 1; tol = 1e-2;

while flag

    if iS > 2
        clf( fig_ri )
        [~,ax_siz,Leg_siz,idx_chs] = sizing_plot(Des_Air.TLARs,CLmax_TO_vett,...
            CLmax_LND_vett,CLmax_CR_vett,Des_Air.MLndoMTo,sigma,...
            Des_Air.SizHis(iS-1).CD0,Des_Air.dCD0_wave,TisaoT50,Des_Air.MCroMTo,...
            V_cr_vet,h_cr_vet,phi_v,fig_ri,fig_aux,ch_idxs);
    end
    figure(fig_ri.Number);
    %subplot 211
    % Plots a line corresponding to the assumed WoS
    fig_it           = plot(Des_Air.SizHis(iS-1).WoS*[1,1],[0,1]);
    fig_it.LineStyle = '--'; fig_it.LineWidth = 2; %fig(1,i).Marker = aero_obj(i).Mark; fig(1,i).MarkerSize = 4;
    col              = rand(1,3); fig_it.MarkerEdgeColor = col;
    %axes
    %legend( fig_it,['(W/S)_{it ',num2str(iS-1),'} = ',num2str(Sizing(iS-1).WoS)],'Interpreter','Latex' );
    %legend( [ax_siz,fig_it],{Leg_siz{1:end},['(W/S)_{it ',num2str(iS-1),'} = ',num2str(Sizing(iS-1).WoS)]},'Interpreter','Latex' );
    
    tmp = input('Choose W/T');
    if tmp == -1 && iS > 2
        Des_Air.SizHis(iS-1).ToW = Des_Air.SizHis(iS-2).ToW;
    else
        Des_Air.SizHis(iS-1).ToW = tmp;
    end
    figure(fig_ri.Number)
    fig_pt = plot( Des_Air.SizHis(iS-1).WoS,Des_Air.SizHis(iS-1).ToW );
    fig_pt.LineStyle = 'none';  fig_pt.Marker = 'o'; fig_pt.MarkerSize = 6;
    col = rand(1,3); fig_pt.MarkerEdgeColor = col;
    %legend( [ax_siz,fig_it,fig_pt],{Leg_siz,['(W/S)_{it ',num2str(iS-1),'} = ',num2str(Sizing(iS).WoS)],'Sizing Point'},'Interpreter','Latex' );

    tmp = input('Choose W/S');
    if tmp == -1
        Des_Air.SizHis(iS).WoS = Des_Air.SizHis(iS-1).WoS;
    else
        Des_Air.SizHis(iS).WoS = tmp;
    end
    Des_Air.SizHis(iS).S                    = Des_Air.MTOM_est0/Des_Air.SizHis(iS).WoS;
    [Des_Air.SizHis(iS).CD0,Des_Air.SizHis(iS).Swet] = Des_Air.polar_est( Des_Air.SizHis(iS).S,Des_Air.MTOM_est0 ); %polar_est(Sizing(iS).S,MTOM_it0);
    
    err.WoS = abs( (Des_Air.SizHis(iS).WoS - Des_Air.SizHis(iS-1).WoS)/Des_Air.SizHis(iS-1).WoS ) ;
    %abs.WoT = ( Sizing(iS).WoT - Sizing(iS-1).WoT ) ;
    flag = err.WoS > tol; %&& abs.WoT <tol;
    iS = iS + 1;
    %hold off

end

sizing_plot_cfr(Airl,nAero,fig_ri,ax_siz,Leg_siz,fig_pt)
%

%% Wing Design
Des_Air = Des_Air.final_out( idx_chs,CLmax_TO_vett,CLmax_CR_vett,...
    CLmax_LND_vett,V_cr_vet,h_cr_vet );
Wing_Design( Des_Air );
