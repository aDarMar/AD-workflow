function Mdd = buffet_check(M_vett,Cl_max,tc_mean,cosc4,deltamcc)
% Mdd Fa il check sulla buffet barrier
% 
% INPUT:
%   M_vett   - Vettore di mach assunto 
%   Cl_max   - Cl max in clean condition (M=0)
%   tc_mean  - t/ medio 
%   cosc4    - cos dell'angolo di freccia a c/4
%   deltamcc - delta mach critico per profili supercritici(0.06)
%
% OUTPUT:
%   
for i=1:length(M_vett)
    Cl_d0(i) = Cl_max/sqrt(1-M_vett(i)^2);% CL con correzione di prandtl-glauert
end
x  = (tc_mean/cosc4));
K1 = 2.8355*x^2-1.9072*x+0.9499;
K2 = 0.2*(1-2.131*x);
for i=1:length(M_vett)
    Cl_mcc = (K1*cosc4^2-(M_vett(i)-deltamcc)*cosc4^3)/K2;
end
diff = Cl_mcc-Cl_d0;
jstart = find(diff<0.01,1,"first");
Msub = M_vett(jstart:end);
for j = 1:length(Msub)
            Mdd(j) = (Msub(j)-0.06)/(1.02+0.08*(1-cosc4));
            Cl_mdd(j) = (K1*cosc4^2-Mdd(j)*cosc4^3)/K2;
end
mc = 0.785; cl_c= 0.38;%current point
plot(M_vett(1:jstart),Cl_d0(1:jstart));hold on;
plot(Msub,Cl_mcc(jstart:end));hold on;
plot(Msub,Cl_mdd);hold on;
plot(mc,cl_c,'o');
legend('Cl a M diverso da 0','Cl a Mcc','Cl a MDD','Current point')

