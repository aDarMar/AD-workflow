function [MTOM_it0,EM_it0] = Class0_MTOM( air_des,a_coeffs,fig,x0 )
%Class0_MTOM Function that evaluates the intersection between statistical 
% and Linear MTOMvsME curves
%   Detailed explanation goes here

a = a_coeffs(1); b = a_coeffs(2);
ax_3 = subplot(1,3,3,'Parent',fig);
Wempty_reg  = @(Wmtom_reg) ( 10.^( ( log10(Wmtom_reg*2.2046) - a )./b )/2.2046 ) ;      % in [kg]
Wempty_stat = @(Wmtom_reg) air_des.weight_coeff_c*Wmtom_reg-air_des.weight_coeff_d;
find_W      = @(Wmtom_reg) Wempty_reg(Wmtom_reg) - Wempty_stat(Wmtom_reg);
if nargin <4
    x0          = 120000; % [Kg]
end
MTOM_it0    = fzero( find_W,x0 );
EM_it0      = Wempty_reg( MTOM_it0 );

Wmtom_reg   = (2:20)*1e4;
a1 = plot( ax_3,Wmtom_reg,Wempty_reg( Wmtom_reg ) ); hold( ax_3,'on' );
a2 = plot( ax_3,Wmtom_reg,Wempty_stat( Wmtom_reg ) );
ax_3.XGrid = 'on'; ax_3.YGrid = 'on'; %axis(ax_2,'equal');
ax_3.XMinorGrid = 'on'; ax_3.YMinorGrid = 'on'; 


end

