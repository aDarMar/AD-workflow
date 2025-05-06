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

% for i=1:rows
%     for k=1:column
%         if k==1
%             A(i,k)=1;
%         else
%             A(i,k)=log10(WE(i));
%         end %close if
%     end     % close i for
% end         % close k for
% 
% % Known vector b
% b=log10(WTO)';
% 
% % Pseudoinverse matrix calculation and solution of the system
%  solution = pinv(A)*b;
%  a=solution(1);
%  b=solution(2);



end

