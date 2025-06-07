function [ ax1,ax2 ] = sizing_plot_TO( ax1,ax2,Smax_TO,CLmaxes,sigma )
%sizing_plot_TO Function that draws the limitation on Take-Off for a jet
%airplane
%INPUT
%   ax1: axis object to the W/S in [kg/m^2]
%   ax1: axis object to the W/S in [lb/ft^2]
%   Smax_TO: Take-Off field Length (FAR25) in [m]
%   CLmaxes: vector of selected CLmaxes in Take-Off
%   sigma : rho/rhoSL for the given airport

colors = [
    0.89, 0.10, 0.11;   % rosso acceso
    1.00, 0.50, 0.00;   % rosso aranciato
    0.70, 0.13, 0.13;   % rosso mattone
    0.98, 0.50, 0.45;   % rosso rosa
    0.55, 0.00, 0.00;   % rosso scuro
    0.82, 0.00, 0.00    % rosso ciliegia
];

WoS = linspace(0,1000,2);%*0.4535924/(0.3048^2);    %[Kg/m^2]
nCLmaxes = length( CLmaxes );
if nCLmaxes > length( colors )
    aux_c = random( nCLmaxes-colors,3 );
    colors = [colors;aux_c];
end
i = 1;
ToW = @(WoS,CLmax) 37.5/( sigma*CLmax*Smax_TO*3.28084 )*( WoS*2.204623/(3.28084^2) );
% Kg/M^2
lin(i,1) = plot( ax1,WoS,ToW(WoS,CLmaxes(i)) ); hold( ax1,'on'); % !!!! 
% hold on hold on si applica solo all’axes corrente, cioè a gca.
lin(i,1).LineStyle = '--'; lin(i,1).LineWidth = 1.5; lin(i,1).Color = colors(i,:);
% lb/ft^2
lin(i,2) = plot( ax2,WoS*2.204623/(3.28084^2),ToW(WoS,CLmaxes(i)) ); hold( ax2,'on');
lin(i,2).LineStyle = '--'; lin(i,2).LineWidth = 1.5; lin(i,2).Color = colors(i,:);


LEG{1} = ['Take-Off with CL$_{max,TO}$ = ', num2str(CLmaxes(i))];
if nCLmaxes > 1
    for i = 2:nCLmaxes
        %subplot 211
        lin(i,1) = plot( ax1,WoS,ToW(WoS,CLmaxes(i)) );
        lin(i,1).LineStyle = '--'; lin(i,1).LineWidth = 1.5; lin(i,1).Color = colors(i,:);
        %set(fig(i,1), 'Visible','off');
        %subplot 212
        lin(i,2) = plot( ax2,WoS*2.204623/(3.28084^2),ToW(WoS,CLmaxes(i)) );
        lin(i,2).LineStyle = '--'; lin(i,2).LineWidth = 1.5; lin(i,2).Color = colors(i,:);
        %set(fig(i,2), 'Visible','off');
        LEG{i} = ['Take-Off with CL$_{max,TO}$ = ', num2str(CLmaxes(i))];
    end
end
%subplot 211
legend( ax1,LEG,'Interpreter','latex','FontSize',16 );
%subplot 212
legend( ax2,LEG,'Interpreter','latex','FontSize',16 );
end