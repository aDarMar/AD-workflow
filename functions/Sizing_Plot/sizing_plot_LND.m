function [ ax1,ax2 ] = sizing_plot_LND( ax1,ax2,Smax_LND,CLmaxes,sigma,WLND2WTO )
%sizing_plot_LND Function that draws the Landing limitation for a jet
%airplane
%   Smax_LND: design (not FAR25) length  in [m]
%   CLmaxes: vector of selected CLmaxes in Landing
%   sigma : rho/rhoSL for the given airport
%   WLND2WTO: Landing weight as a fraction of Take-Off Weight
nCLmaxes = length( CLmaxes );

%% Graphics
colors = [
    0.106, 0.620, 0.467;  % Verde brillante
    0.196, 0.804, 0.196;  % Verde lime
    0.333, 0.420, 0.184;  % Verde oliva scuro
    0.596, 0.984, 0.596;  % Verde menta
    0.133, 0.545, 0.133;  % Verde foresta
    0.000, 0.392, 0.000;  % Verde scuro intenso
    ];

if nCLmaxes > length( colors )
    aux_c = random( nCLmaxes-colors,3 );
    colors = [colors;aux_c];
end

%% Plots
%cf = 0.3/3.281 * (1/0.5144)^2;
i = 1;

%WoS = @(CLmax) WLND2WTO*3.281^2/2.205*0.5*1.225*sigma*( Smax_LND/(cf*1.3^2) )*CLmax;
% [kg/m^2]
lin(i,1) = plot( ax1,[1,1]*WoS_fun(Smax_LND,CLmaxes(i),sigma,WLND2WTO),[0,1] ); hold( ax1, 'on')
lin(i,1).LineStyle = ':'; lin(i,1).LineWidth = 1.5; lin(i,1).Color = colors(i,:);
lin(i,1).DisplayName = ['Landing with CL$_{max,LND}$ = ', num2str(CLmaxes(i))];

%set(fig(i,1), 'Visible','off');
LEG{1} = ['Landing with $CL_{max,LND}$ = ', num2str(CLmaxes(i))];
% [lb/ft^2]
lin(i,2) = plot( ax2,[1,1]*WoS_fun(Smax_LND,CLmaxes(i),sigma,WLND2WTO)*2.204623/(3.28084^2),[0,1] ); hold( ax2, 'on' );
lin(i,2).LineStyle = ':'; lin(i,2).LineWidth = 1.5; lin(i,2).Color = colors(i,:);
lin(i,2).DisplayName = ['Landing with CL$_{max,LND}$ = ', num2str(CLmaxes(i))];

%set(fig(i,2), 'Visible','off');

if nCLmaxes > 1
    for i = 2:nCLmaxes

        lin(i,1) = plot( ax1,WoS_fun(Smax_LND,CLmaxes(i),sigma,WLND2WTO)*[1,1],[0,1] );
        lin(i,1).LineStyle = ':'; lin(i,1).LineWidth = 1.5; lin(i,1).Color = colors(i,:);
        lin(i,1).DisplayName = ['Landing with CL$_{max,LND}$ = ', num2str(CLmaxes(i))];

        lin(i,2) = plot( ax2,[1,1]*WoS_fun(Smax_LND,CLmaxes(i),sigma,WLND2WTO)*2.204623/(3.28084^2),[0,1] );
        lin(i,2).LineStyle = ':'; lin(i,2).LineWidth = 1.5; lin(i,2).Color = colors(i,:);
        lin(i,2).DisplayName = ['Landing with CL$_{max,LND}$ = ', num2str(CLmaxes(i))];

    end
end
%legend(
%ax1.Legend.String       = [ax1.Legend.String,LEG]; % Updates the Legend string names
%ax1.Legend.PlotChildren = [ax1.Legend.PlotChildren;lin(1,:)'];

%ax2.Legend.String = [ax1.Legend.String,LEG];

%);%,'Interpreter','latex','FontSize',16 );
%legend( ax2,[ax1.Legend.String,LEG] );%,'Interpreter','latex','FontSize',16 );
end

function WoS = WoS_fun(SLand,CLmax,sigma,WLND2WTO)
% WoS in Kg/m^2
Va_sq = (SLand/0.6)*3.281/0.3; % Approach Speed [kts^2]
Va_sq = Va_sq/1.3^2;    % Vstall [kts^2]
Va_sq = Va_sq/(1.944)^2; %Vstall [m/s]^2
WoS = 0.5*1.225*sigma*Va_sq*CLmax/WLND2WTO; %[N/m^2]
WoS = WoS/9.81; %[kg/m^2]
%WoS*2.205/(3.281^2);
%WLND2WTO*3.281^2/2.205*0.5*1.225*sigma*( Va_sq/1.3^2 )*CLmax;
end
