function [ ax1,ax2 ] = sizing_plot_Cruise( ax1,ax2,CD0,dCD0_wave,Vcr,h_cruise,n_fact,WcroWTO,e,AR,phi_v,T0oTcr)
%sizing_plot_Cruise Function that draws the cruise limitation for a jet
%airplane. Inputs can be vector and in such case, every element corresponds
%to a specific flight condition
%   CD0: zero-lift drag coefficie in cruise (clean config)
%   dCD0_wave: wave drag addition
%   Vcr: cruise speed in [m/s]
%   h_cruise : cruise altitude in [m]
%   WTOoWcr: vector containing weigth ratios between [end of climb (4), end
%       of cruise (5) ] weights over WTO weight
%   e: oswald factor in cruise Uclean config.)
%   AR: Imposed AR
%   phi_v: vector of admission rates [-]
%   T0oTcr: ratio of max static thrust over Thrust in cruise. If not given
%       the model Tcr = T0*phi*sigma will be used

nConds = length( Vcr );
%% Graphics
colors = [
    0.22, 0.49, 0.72;  % blu brillante
    0.53, 0.81, 0.92;  % azzurro chiaro
    %0.27, 0.51, 0.71;  % blu acciaio
    0.00, 0.74, 0.83;  % azzurro-verde (cyan)
    0.00, 0.00, 0.55;  % blu scuro
    0.42, 0.35, 0.80;  % blu-grigio
    ];

if nConds > length( colors )
    aux_c = random( nCLmaxes-colors,3 );
    colors = [colors;aux_c];
end
MARK = {'x','^','square'};
%% Plot
COND = {'Initial W ','Final W ','Avg. W '};
nW = length( WcroWTO );
if nW == 2
    WcroWTO(3) = ( WcroWTO(1) + WcroWTO(2) ) *0.5;
    nW = nW + 1;
end
K = 1/(pi*AR*e);
WoS = linspace(0,1000,100); % [Kg/m^2]
for j = 1:nConds
    [T, a, P, rho] = atmosisa(h_cruise(j));
    
    sigma = rho/1.225;
    if nargin < 12
        T0oTcr = 1/(0.71*sigma*phi_v(j));
    end
    q = 0.5*sigma*1.225*Vcr(j)^2;

    for i = 1:nW
        ToW = ( (CD0+dCD0_wave)*q./(WoS*9.81)/WcroWTO(i) + K/q .* (WoS*9.81)*WcroWTO(i)*n_fact(j)^2 )*WcroWTO(i)*T0oTcr; %WoS in [Kgf]
        % [kg/m^2]
        lin(2*i-1,j) = plot( ax1,WoS,ToW );
        lin(2*i-1,j).LineStyle = '-'; lin(2*i-1,j).LineWidth = 2;
        lin(2*i-1,j).Color = colors(j,:); lin(2*i-1,j).Marker = MARK{i};
        lin(2*i-1,j).MarkerSize = 2.5;
        lin(2*i-1,j).DisplayName = [COND{i},[' at h = ',num2str( h_cruise(j) ),...
            ' V = ',num2str( Vcr(j) ),' $\phi$ = ',num2str( phi_v(j) ),' n = ',num2str(n_fact(j))] ];

        % [lb/ft^2]
        lin(2*i,j) = plot( ax2,WoS*2.204623/(3.28084^2),ToW );
        lin(2*i,j).LineStyle = '-'; lin(2*i,j).LineWidth = 2;
        lin(2*i,j).Color = colors(j,:); lin(2*i,j).Marker = MARK{i};
        lin(2*i,j).MarkerSize = 2.5;
        lin(2*i,j).DisplayName = [COND{i},[' at h = ',num2str( h_cruise(j) ),...
            ' V = ',num2str( Vcr(j) ),' $\phi$ = ',num2str( phi_v(j) ),' n = ',num2str(n_fact(j))] ];
    end
    
end
end
