function [fig,LEG] = sizing_plot_TO( Smax_TO,CLmaxes,sigma )
%sizing_plot_TO Function that draws the limitation on Take-Off for a jet
%airplane
%   Smax_TO: Take-Off field Length (FAR25) in [m]
%   CLmaxes: vector of selected CLmaxes in Take-Off
%   sigma : rho/rhoSL for the given airport

WoS = linspace(0,1000,2);%*0.4535924/(0.3048^2);    %[Kg/m^2]
nCLmaxes = length( CLmaxes );


i = 1;
ToW = @(WoS,CLmax) 37.5/( sigma*CLmax*Smax_TO*3.28084 )*( WoS*2.204623/(3.28084^2) );
subplot 211
fig(i,1) = plot( WoS,ToW(WoS,CLmaxes(i)) ); hold on
LEG{1} = ['Take-Off with $CL_{max,TO}$ = ', num2str(CLmaxes(i))];
fig(i,1).LineStyle = '--'; fig(i,1).LineWidth = 1.5;
%set(fig(i,1), 'Visible','off');
subplot 212
fig(i,2) = plot( WoS*2.204623/(3.28084^2),ToW(WoS,CLmaxes(i)) ); hold on
fig(i,2).LineStyle = '--'; fig(i,2).LineWidth = 1.5;
%set(fig(i,2), 'Visible','off');
if nCLmaxes>1
    for i = 2:nCLmaxes
        subplot 211
        fig(i,1) = plot( WoS,ToW(WoS,CLmaxes(i)) );
        fig(i,1).LineStyle = '--'; fig(i,1).LineWidth = 1.5;
        %set(fig(i,1), 'Visible','off');
        subplot 212
        fig(i,2) = plot( WoS*2.204623/(3.28084^2),ToW(WoS,CLmaxes(i)) );
        fig(i,2).LineStyle = '--'; fig(i,2).LineWidth = 1.5;
        %set(fig(i,2), 'Visible','off');
        LEG{i} = ['Take-Off with CL$_{max,TO}$ = ', num2str(CLmaxes(i))];
    end
end
subplot 211
legend( fig(:,1),LEG,'Interpreter','latex','FontSize',16 );
subplot 212
legend( fig(:,2),LEG,'Interpreter','latex','FontSize',16 );
end