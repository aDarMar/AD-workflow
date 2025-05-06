function [Cd0,Swet] = polar_est(S,MTOM)
%polar_est Function that evaluates the polar with a statistical approach.
%   MTOM:   in Kg
%   S   :   estimated wing surface [m]
global main_fold

lb2kg = 0.45359237; ft2m = 0.3048;
% From Roskam, values for Transport Jets
c = 0.0199; d = 0.7331;

Swet  = 10^( c+d*log10(MTOM/lb2kg) ); %Swet in ft^2
temp = readmatrix([main_fold,'\statistical_data\cf_vs_Swet.csv']);
Cf_eq = @(Sw) interp1(temp(:,1),temp(:,2),Sw);
a   = [ 2.0458,2.0969,2.1549,2.2218,2.301,2.3979,2.5229,2.699 ]*(-1); a = flip(a);
b =  ones( 1,length(a) );
Cfs = (2:9)*1e-3;

a_f =@(Cf) a(1)*(Cf<Cfs(1)) + interp1(Cfs,a,Cf,'linear',0) ... %*( ~( ~(Cf<Cfs(1) )&& ~( Cf>Cfs(end) ) ) )...
    + a(end)*(Cf>Cfs(end));
b_f =@(Cf) b(1)*(Cf<Cfs(1)) +interp1(Cfs,b,Cf,'linear',0) + b(end)*(Cf>Cfs(end));

f = 10^(a_f(Cf_eq(Swet) ) + b_f(Cf_eq(Swet))*log10( Swet )); % f in [ft^2]
f = f*(ft2m^2); % f in [m^2]

Cd0 = f/S;

Swet = Swet/(ft2m^2);

end


