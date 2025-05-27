function [MTOM_it0,EM_it0] = Class0_MTOM( air_des,a_coeffs,x0 )
%Class0_MTOM Function that evaluates the intersection between statistical 
% and Linear MTOMvsME curves
%   Detailed explanation goes here

a = a_coeffs(1); b = a_coeffs(2);

Wempty_reg  = @(Wmtom_reg) ( 10.^( ( log10(Wmtom_reg*2.2046) - a )./b )/2.2046 ) ;      % in [kg]
Wempty_stat = @(Wmtom_reg) air_des.weight_coeff_c*Wmtom_reg-air_des.weight_coeff_d;
find_W      = @(Wmtom_reg) Wempty_reg(Wmtom_reg) - Wempty_stat(Wmtom_reg);
if nargin <3
    x0          = 120000; % [Kg]
end
MTOM_it0    = fzero( find_W,x0 );
EM_it0      = Wempty_reg( MTOM_it0 );

figure(); hold on
Wmtom_reg   = (2:20)*1e4;
a1 = plot( Wmtom_reg,Wempty_reg( Wmtom_reg ) );
a2 = plot( Wmtom_reg,Wempty_stat( Wmtom_reg ) );
grid minor

end

