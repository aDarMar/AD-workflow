function L = L_funct(b,c,sweep,y,s_y)
%L_funct: Weissinger L-Function from the same method. Formulation taken from
%NACA Report No. 921
%   b: wing span
%   c: chord at eta section
%   sweep: sweep angle at c/4 in deg
%   y: spanwise coordinate at which the chord is referred
%   s_y: spanwise coordinate of the vortex filament 
%OUTPUT
%   L: Weissinger L-function value for y,y_s. From the reference it is: 
%   ---------------------------- L(y,s_y) ---------------------------- 
boc = b/c;  tS4 = tan( sweep*pi/180 );
eta = 2*y/b; s_eta = 2*s_y/b; a_eta = abs(eta);
if eta < 0
    L = 1/( boc*(eta-s_eta) )*( sqrt( ( 1+boc*(a_eta+s_eta)*tS4 )^2 + boc^2*(eta-s_eta)^2 )/( 1+boc*(a_eta+eta)*tS4 ) - 1 ) + ...
        ( 2*tS4*sqrt( (1+boc*a_eta*tS4)^2 + (boc*eta)^2 ) )/( (1+boc*(a_eta-eta)*tS4)*(1+boc*(a_eta+eta)*tS4) );
else
    L = 1/( boc*(eta-s_eta) )*( sqrt( ( 1+boc*(a_eta-s_eta)*tS4 )^2 + boc^2*(eta-s_eta)^2 )/( 1+boc*(a_eta-eta)*tS4 ) - 1 );
end

end

