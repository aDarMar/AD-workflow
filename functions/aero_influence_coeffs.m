function a_coeffs = aero_influence_coeffs( m,M,geom_vec,aero_vec,b )
%aerosymload: function finds the circulation solving the circulation at m
%spanpoints
%   c: chords given at exactly Multhopp integration points
%   M: number of points chosen to perform the numerical integration of L

m = length( c );
[phi_v,eta]     = phi_funct(m);     %vector of angles and non-dimensional positions of stations along the span
[phi_int,b_eta] = phi_funct(M);     %vector of angles and non-dimensional positions of integration points along the span
% Building a section class
des_wing = PaneledWing( m,M,geom_vec,aero_vec,b );
m_red           = (m+1)/2;          %Number of reduced points
a_coeffs        = nan(m_red);       %Matrix of symmetric influence coefficients

for nu = 1:mred
    % Builds the matrix of coefficients row by row
    a_coeffs(nu,:) = aerosymmbuilder( nu,m,M,phi_int,b_eta,phi_v,eta,c_vect(nu),b );
end
    
    
end

function a_nu = aerosymmbuilder(nu,m,M,phi_int,b_eta,phi_vect,eta,c,b)
%aerosymmbuilder: function that builds the equation for circulation in case of
%symmetric load. In other words, it builds the nu-th row of the influence
%coefficients matrix.
%   nu: nu-th row of equation. In physical terms it represents the index of
%       the nu-th control point in which we are imposing the no flow
%       through condition. NU IS A SCALAR
%   m: number of such control points
%   M: number of points of integrations chosen for the integration of L
%   phi_int, b_eta: angular and non dimensional coordinates of the
%       aforementioned points
%   phi_vect,eta: angular and non dimensional coordinates of the spanwise
%       stations
%   c: choord of nu (control point) section. c_nu 
%   b: wing span

m_red   = 0.5*(m+1);
a_nu    = nan( 1,m_end );
n_idxs  = 1:m_red-1; 
n_idxs  = n_idxs( n_idxs~=nu );

% For-cycle from 1 to (m+1)/2 - 1 excluding nu
for n = n_idxs
    % Eq. (A37) case n =/= nu
    B       = littlebfun( nu,n,phi_vect ) + littlebfun( nu,m+1-n,phi_vect );
    a_nu(n) = -2*B + b/c*gbarfun( M,m,nu,phi_int,b_eta,phi_vect(n),eta,b,c );
end
% Adding (m+1)/2 point: this step is made outside the for cycle because the
% functions gbar and B assume  different values for n = (m+1)/2
B          = littlebfun( nu,i,phi_vect );
a_nu(i)    = -2*B + b/c*gbarfun();
% Adding nu point
a_nu(nu)   = 2*(m+1)/(4*sin(phi_vect(nu)));    %b_nu,nu expression
if nu == m_red
     a_nu(nu) = a_nu(nu) + gbarfun_special( M,m,nu,phi_int,b_eta,phi_vect,eta,b,c );
else
    % The function g
    a_nu(nu) = a_nu(nu) + gbarfun( M,m,nu,phi_int,b_eta,phi_vect(nu),eta,b,c );
end


end


function b_ij = littlebfun( nu,n,phi_vect )
%littlebfun: function that evaluates the coefficient b(nu,n)
%   nu: index of the control point
%   phi_n: index of the other point
%   phi_vect: vector of angular positions of chosen sections
    m = length( phi_vect );
    b_ij = sin( phi_vect(n) )/( ( cos(phi_vect(n))-cos(phi_vect(nu)) )^2 )*( (1-(-1)^(n-nu))/(2*(m+1)) );
end

function [phi_n,eta] = phi_funct( m )
%phi_funct: function that calculates spanwise sections according to
%Multhopp quadrature method
%   m: number of points
phi_n = nan(m,1); eta = phi_n;
for i = 1:m
    phi_n(i) = n*pi/(m+1);
    eta(i) = cos(phi_n(i));
end
end

function gval = gbarfun( M,m,nu,phi_int,b_eta,phi_n,eta,b,c,Sweep )
%gbarfun: function that evaluates g_bar in case the second index is
%different from (m+1)/2. 
% ------------------ BE CAREFUL! ------------------
%In case the value of g for that index is required,
%refer to the function gbarfun_special
%--------------------------------------------------
%   M: number of integration points for L
%   phi_n: angular coordinate of the spanwise station
%   phi_int: vector of angular coordinates of integration points for L
%   m: number of spanwise coordinates

% Case mu = 0 and n =/= (m+1)/2
gval = ffun( m,phi_n,0 )*L_funsymm( eta(nu),1,b,c,Sweep );
Mred = (M-1)/2;
for mu = 1:Mred
    % Case mu =/= 0 and n =/= (m+1)/2
    gval = gval + 2*ffun( m,phi_n,phi_int(mu) )*L_funsymm( eta(nu),b_eta(mu),b,c,Sweep );
end
gval = gval*(-1)/(2*(M+1));
end

function gval = gbarfun_special( M,m,nu,phi_int,b_eta,phi_n,eta,b,c,Sweep )
%gbarfun: function that evaluates g_bar in case the second index IS
%EQUAL TO (m+1)/2. 
% ------------------ BE CAREFUL! ------------------
%In case the value of g for the other indexes is required,
%refer to the function gbarfun
%--------------------------------------------------
%   M: number of integration points for L
%   phi_n: angular coordinate of the spanwise station
%   phi_int: vector of angular coordinates of integration points for L
%   b_eta: vector of non-dimensional coordinates of integration points
%   m: number of spanwise coordinates

% Case mu = 0 and n == (m+1)/2
gval = 0.5*ffun( m,phi_n,0 )*L_funsymm( eta(nu),1,b,c,Sweep );
Mred = (M-1)/2;
for mu = 1:Mred
    % Case mu =/= 0 and n == (m+1)/2
    gval = gval + ffun( m,phi_n,phi_int(mu) )*L_funsymm( eta(nu),b_eta(mu),b,c,Sweep );
end
gval = gval*(-1)/(2*(M+1));
end

function fnm = ffun( m,phi_n,phi_mu )
    for mu1 = 1:m
        fnm = fnm + mu1*sin( mu1*phi_n )*cos( mu1*phi_mu );
    end
    fnm = fnm*2/(m+1);
end

function L = L_funsymm( eta,b_eta,b,c) % CONTROLLA LE CHIAMATEEEEEEE
    boc = b/c; 
    L = 1/( boc*(eta-b_eta) )*( sqrt( (1+boc*(eta-b_eta)*tS4)^2 + (boc*(eta-b_eta)^2) ) - 1 ) -...
        1/( boc*(eta-b_eta) )*( sqrt( (1+boc*(eta-b_eta)*tS4)^2 + (boc*(eta-b_eta)^2) )/( 1+2*boc*eta*tS4 ) - 1 ) -...
        ( 2*tS4*sqrt( (1+boc*eta*tS4)^2+(boc*eta)^2 ) )/(1+2*boc*eta*tS4);
end