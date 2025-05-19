function a_coeffs = aero_influence_coeffs( m,M,geom_vec,aero_vec,b,sweep )
%aerosymload: function finds the circulation solving the circulation at m
%spanpoints
%   c: chords given at exactly Multhopp integration points
%   M: number of points chosen to perform the numerical integration of L

%m = length( geom_vec(:,1) );
% Building a section class
des_wing = PaneledWing( m,M,geom_vec,aero_vec,b,sweep );
m_red           = (m+1)/2;          %Number of reduced points
a_coeffs        = nan(m_red);       %Matrix of symmetric influence coefficients

for nu = 1:m_red
    % Builds the matrix of coefficients row by row
    a_coeffs(nu,:) = aerosymmbuilder( nu,des_wing );
end

end

function a_nu = aerosymmbuilder( nu,dwing )
%aerosymmbuilder: function that builds the equation for circulation in case of
%symmetric load. In other words, it builds the nu-th row of the influence
%coefficients matrix.
%   nu: nu-th row of equation. In physical terms it represents the index of
%       the nu-th control point in which we are imposing the no flow
%       through condition. NU IS A SCALAR
%   des_wing: a PaneledWing object containing all geometrical data for the
%       wing.

m_red   = 0.5*(dwing.m+1);
a_nu    = nan( 1,m_red );
n_idxs  = 1:m_red-1; 
n_idxs  = n_idxs( n_idxs~=nu );

% For-cycle from 1 to (m+1)/2 - 1 excluding nu
for n = n_idxs
    % Eq. (A37) case n =/= nu
    B       = dwing.littlebfun( nu,n ) + dwing.littlebfun( nu,dwing.m+1-n );
    a_nu(n) = -2*B + dwing.b/dwing.geom_sect(nu).c*dwing.gbarfun( nu,n );
end
% Adding (m+1)/2 point: this step is made outside the for cycle because the
% functions gbar and B assume  different values for n = (m+1)/2
B          = dwing.littlebfun( nu,m_red );
a_nu(m_red)= -2*B + dwing.b/dwing.geom_sect(nu).c*dwing.gbarfun_special( nu,m_red );
% Adding nu point
a_nu(nu)   = 2*(dwing.m+1)/( 4*sin(dwing.geom_sect(nu).phi) );    %2*b_nu,nu expression
if nu == m_red
    % case nu = n = (m+1)/2
     a_nu(nu) = a_nu(nu) + dwing.b/dwing.geom_sect(nu).c*dwing.gbarfun_special( nu,nu );
else
    % case nu = n =/= (m+1)/2
    a_nu(nu) = a_nu(nu) + dwing.b/dwing.geom_sect(nu).c*dwing.gbarfun( nu,nu );
end


end
