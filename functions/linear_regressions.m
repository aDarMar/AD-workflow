function [a_coeff,b_coeff] = linear_regressions(aero_obj,nAero)
%linear_regressions Function that evaluates a,b coefficient of the linear
%regression between WTO and W_Empty.
%   log10(WTO) = a + b*log10(W_empty)
%   with WTO, W_empty in [lb]
%   OUTPUT
%   a_coeff
%   b_coeff
% Unknowns matrix A
lb2kg = 0.45359237;

A = ones(nAero,2); b = nan(nAero,1);
for i = 1:nAero
    % Masse in libbre
    A(i,2) = log10( aero_obj(i).EM/lb2kg ); 
    b(i)   = log10( aero_obj(i).MTOM/lb2kg );
end

sol = pinv(A)*b;
a_coeff = sol(1) ; b_coeff = sol(2);


end

