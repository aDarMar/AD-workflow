close all; clear; clc


M = 7; m_in = 2; m = 7;
geom_vec = nan(m_in,9); aero_vec = nan(m_in,9);
aero_vec(:,1) = 0; aero_vec(:,2) = 2*pi;
%temp = [ 3.152240935; 3.585786438; 4.234633135;5 ];
%geom_vec(:,1) = [ temp;flip(temp(1:end-1)) ];
geom_vec(:,1) = [1,0]; geom_vec(:,2) = [3,5];

b = 32; sweep = 20;
aero_influence_coeffs( m,M,geom_vec,aero_vec,b,sweep )