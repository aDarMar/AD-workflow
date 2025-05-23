function [outputArg1,outputArg2] = basicLoad(dwing,inputArg2)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here


% Reduced Matrix for Basic Loading
% i: i-th control point (nu in reference)
% j: j-th element
% r: index of reference point. IT IS ASSUMED THIS IS THE ROOT, otherwise
% the equation below is no longer valid
% A(i,j) = a(i,j)-a(r,j)-( a(i,r)-a(r,r) )*2*sin( phi(j) )
nu_r = (dwing.m+1)*0.5; % reference station is the root chord
nu_idx = 1:(dwing.m+1)*0.5; nu_idx_red = nu_idx( nu_idx ~= nu_r);
A = nan( (dwing.m+1)*0.5 - 1 );
for nu = nu_idx_red
    for n = nu_idx
        A(nu,n) = dwing.a_coeffs(nu,n)-a_coeffs(nu_r,n)-...
            ( a_coeffs(nu,nu_r)-a_coeffs(nu_r,nu_r) )*2*sin( dwing.geom_sect(nu).phi );
    end
end
Gs = A\eps;
if nu_r > 1
    Gs = [Gs(1:nu_r-1);0,Gs(nu_r:end)];
else
    Gs = [0;Gs];
end
a0L_ref = 0;
for n =1:nu_idx_red
    % Equation for total CL = 0
    Gs(nu_r) = Gs(nu_r) - Gs(n)*2*sin( dwing.geom_sect(nu).phi );
    % Equation for circulation at Reference Station
    a0L_ref  = a0L_ref + ( a_coeffs(nu_r,n) - a_coeffs(nu_r,nu_r)...
        *2*sin( dwing.geom_sect(n).phi ) )*Gs(n);
end

end