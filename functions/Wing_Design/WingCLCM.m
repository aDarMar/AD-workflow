function [CLw_M0,CLw_Mc,CLw_Ms,CMw] = WingCLCM(yroot,ykink,ytip,croot,ckink,ctip, ...
    xLE_root,xLE_kink,xLE_tip,eps_root,eps_kink,eps_tip,azl_root,azl_kink, ...
    azl_tip,Xac_root,Xac_kink,Xac_tip,CMac_root,CMac_kink,CMac_tip,Cla_root, ...
    Cla_kink,Cla_tip,Cl_max_root,Cl_max_kink,Cl_max_tip, ...
    Cl_star_root,Cl_star_kink,Cl_star_tip,Sw, ...
    k1,k2,k3,Mcruise,sweep_c2,ARw,atry,Ms,Rapp,delta_a,amin_stall)
%INPUT
% yroot        : posizione della root lungo lo span (0)
% ykink        : posizione della kink lungo lo span
% ytip         : posizione della tip  lungo lo span
% croot        : corda del profilo di root
% ckink        : corda del profilo di kink
% ctip         : corda del profilo di tip
% xLE_root     : posizione del LE della root
% xLE_kink     : posizione del LE della kink
% xLE_tip      : posizione del LE della tip
% eps_root     : twist del profilo di root
% eps_kink     : twist del profilo di kink
% eps_tip      : twist del profilo di tip
% azl_root     : alpha_0Lift del profilo di root
% azl_kink     : alpha_0Lift del profilo di kink
% azl_tip      : alpha_0Lift del profilo di tip
% Xac_root     : posizione AC del profilo di root
% Xac_kink     : posizione AC del profilo di kink
% Xac_tip      : posizione AC del profilo di tip
% CMac_root    : Cm del profilo di root rispetto al centro aerodinamico
% CMac_kink    : Cm del profilo di kink rispetto al centro aerodinamico
% CMac_tip     : Cm del profilo di tip rispetto al centro aerodinamico
% Cla_root     : Cl_alpha del profilo di root
% Cla_kink     : Cl_alpha del profilo di kink
% Cla_tip      : Cl_alpha del profilo di tip
% Cl_max_root  : Cl_max del profilo di root
% Cl_max_kink  : Cl_max del profilo di kink
% Cl_max_tip   : Cl_max del profilo di tip
% Cl_star_root : Cl* del profilo di root
% Cl_star_kink : Cl* del profilo di kink
% Cl_star_tip  : Cl* del profilo di tip
% Sw           : superficie alare
% k1           : coeffiente per il triangolo 1 dell'ala
% k2           : coeffiente per il triangolo 2 dell'ala
% k3           : coeffiente per il triangolo 3 dell'ala
% Mcruise      : Mach di crociera
% ARw          : aspect ratio alare
% sweep_c2     : angolo di frecca a metà della corda
% atry         : massimo angolo d'incidenza in crociera
% Ms           : Mach di stallo (0.2)
% Rapp         : rapporto tra i CL_max dell'ala e quello del profilo medio (.8)
% delta_a      : distanza tra alpha* e alpha_stallo
% amin_stall   : valore minimo di alpha per costruire la retta del CL
%OUTPUT
% CLw_M0       : CL dell'ala a M = 0
% CLw_Mc       : CL dell'ala al Mach di crociera
% CLw_Ms       : CL dell'ala al Mach di stallo
% CMw          : CM dell'ala

eta_root = yroot/ytip;   xTE_root = xLE_root + croot;
eta_kink = ykink/ytip;   xTE_kink = xLE_kink + ckink;
eta_tip  = ytip/ytip;    xTE_tip  = xLE_tip + ctip;

c_reg_1    = polyfit([yroot;ykink],[croot;ckink],1);
c_reg_2    = polyfit([ykink;ytip],[ckink;ctip],1);
eps_reg_1  = polyfit([yroot;ykink],[eps_root;eps_kink],1);
eps_reg_2  = polyfit([ykink;ytip],[eps_kink;eps_tip],1);
azl_reg_1  = polyfit([yroot;ykink],[azl_root;azl_kink],1);
azl_reg_2  = polyfit([ykink;ytip],[azl_kink;azl_tip],1);
CMac_reg_1 = polyfit([yroot;ykink],[CMac_root;CMac_kink],1);
CMac_reg_2 = polyfit([ykink;ytip],[CMac_kink;CMac_tip],1);
Cla_reg_1  = polyfit([yroot;ykink],[Cla_root;Cla_kink],1);
Cla_reg_2  = polyfit([ykink;ytip],[Cla_kink;Cla_tip],1);
Xac_reg_1  = polyfit([yroot;ykink],[Xac_root;Xac_kink],1);
Xac_reg_2  = polyfit([ykink;ytip],[Xac_kink;Xac_tip],1);
xLE_reg_1  = polyfit([yroot;ykink],[xLE_root;xLE_kink],1);
xLE_reg_2  = polyfit([ykink;ytip],[xLE_kink;xLE_tip],1);

% Tabella Gialla
yvec = 0:1:fix(ytip);
yvec = [yvec,ytip,ykink];       yvec = sort(yvec);
for i=1:length(yvec)
    if yvec(i) < ykink
        cvet(i)     = c_reg_1(2)    + c_reg_1(1)*yvec(i);
        eps_vet(i)  = eps_reg_1(2)  + eps_reg_1(1)*yvec(i);
        azl_vet(i)  = azl_reg_1(2)  + azl_reg_1(1)*yvec(i);
        CMac_vet(i) = CMac_reg_1(2) + CMac_reg_1(1)*yvec(i);
        Cla_vet(i)  = Cla_reg_1(2)  + Cla_reg_1(1)*yvec(i);
        Xac_vet(i)  = Xac_reg_1(2)  + Xac_reg_1(1)*yvec(i);
        xLE_vet(i)  = xLE_reg_1(2)  + xLE_reg_1(1)*yvec(i);
    else
        cvet(i)     = c_reg_2(2)    + c_reg_2(1)*yvec(i);
        eps_vet(i)  = eps_reg_2(2)  + eps_reg_2(1)*yvec(i);
        azl_vet(i)  = azl_reg_2(2)  + azl_reg_2(1)*yvec(i);
        CMac_vet(i) = CMac_reg_2(2) + CMac_reg_2(1)*yvec(i);
        Cla_vet(i)  = Cla_reg_2(2)  + Cla_reg_2(1)*yvec(i);
        Xac_vet(i)  = Xac_reg_2(2)  + Xac_reg_2(1)*yvec(i);
        xLE_vet(i)  = xLE_reg_2(2)  + xLE_reg_2(1)*yvec(i);
    end
    azl_int(i)  = cvet(i)*(azl_vet(i)-eps_vet(i));
    YMAC_int(i) = yvec(i)*cvet(i);
    XMAC_int(i) = xLE_vet(i)*cvet(i);
    xc_4(i)     = Xac_vet(i)*cvet(i)+xLE_vet(i);
    cpower2(i)  = cvet(i)^2;
end

for i=2:length(yvec)
    azl_trap(i) = 0.5*(azl_int(i)+azl_int(i-1))/(yvec(i)-yvec(i-1));
    c2_trap(i)  = 0.5*(cpower2(i)+cpower2(i-1))/(yvec(i)-yvec(i-1));
    YMAC_trap(i)= 0.5*(YMAC_int(i)+YMAC_int(i-1))/(yvec(i)-yvec(i-1));
    XMAC_trap(i)= 0.5*(XMAC_int(i)+XMAC_int(i-1))/(yvec(i)-yvec(i-1));
end

% Quadrato blu
MAC_w     = 2*sum(c2_trap)/Sw;
YMAC_w    = 2*sum(YMAC_trap)/Sw;
XMAC_w    = 2*sum(XMAC_trap)/Sw;
xLE_MAC_w = XMAC_w + xLE_root; 
Xac_w     = xLE_MAC_w + .25*MAC_w;
xTE_MAC_w = xLE_MAC_w + MAC_w;
x1        = Xac_w - xc_4;

% Mean Airfoil
Cla_mean = k1*Cla_tip + k2*Cla_kink + k3*Cla_tip;
Cl_star_mean = k1*Cl_star_root + k2*Cl_star_kink + k3*Cl_star_tip;
Cl_max_mean  = k1*Cl_max_root + k2*Cl_max_kink + k3*Cl_max_tip;
azl_mean = 2*sum(azl_trap)/Sw;
Cla_cruise = (Cla_mean*cosd(sweep_c2))/((sqrt(1-Mcruise^2*(cosd(sweep_c2))^2+( ...
    Cla_mean*cosd(sweep_c2)/pi*ARw)^2))+Cla_mean*cosd(sweep_c2)/pi*ARw);
Cla_Mach = Cla_cruise/sqrt(1-Mcruise^2);

Cm1_int = CMac_vet.*cvet.^2;
Cm2_int = (azl_mean+eps_vet-azl_vet).*Cla_vet.*cvet.*x1;
for i = 2:length(yvec)
    Cm1_trap(i) = 0.5*(Cm1_int(i)+Cm1_int(i-1))*(yvec(i)-yvec(i-1));
    Cm2_trap(i) = 0.5*(Cm2_int(i)+Cm2_int(i-1))*(yvec(i)-yvec(i-1));
end
Cm1 = 2*sum(Cm1_trap)/(Sw*MAC_w);
Cm2 = 2*sum(Cm2_trap)/(Sw*MAC_w);

avet   = -atry:1:atry;  avet = [avet,azl_mean];  avet = sort(avet);
CLw_M0 = Cla_cruise*(avet-azl_mean);
CLw_Mc = Cla_Mach*(avet-azl_mean);
CMw    = Cm1 + Cm2;

% Stallo
Cla_Ms = Cla_cruise/(sqrt(1-Ms^2));
alpha_star = Cla_Ms*azl_mean + Cl_star_mean/Cla_Ms;
CLw_max   = Rapp*Cl_max_mean;
alpha_max = (Cla_Ms*azl_mean+CLw_max)/Cla_Ms+delta_a;

Amat   = [3*alpha_star^2 2*alpha_star 1 0; 3*alpha_max^2 2*alpha_max 1 0;
    alpha_star^3 alpha_star^2 alpha_star 1; alpha_max^3 alpha_max^2 alpha_max 1];
bvet   = [Cla_Ms 0 Cl_star_mean CLw_max]';
sol    = inv(Amat)*bvet;
alpvec = amin_stall:1:alpha_max+1;   alpvec = [alpvec,azl_mean,alpha_star,alpha_max];
alpvec = sort(alpvec);
CLw_Ms = zeros(1,length(alpvec));
for i = 1:length(alpvec)
    if alpvec(i) < alpha_star
        CLw_Ms(i) = Cla_Ms*(alpvec(i)-azl_mean);
    else
        CLw_Ms(i) = sol(1)*alpvec(i)^3 + sol(2)*alpvec(i)^2 + ...
                    sol(3)*alpvec(i) + sol(4);
    end
end
figure()
subplot(1,2,1)
plot(avet,CLw_Mc); xlabel("\alpha"); ylabel("CL"); title(["M=",Mcruise]); grid on
subplot(1,2,2)
plot(alpvec,CLw_Ms); xlabel("\alpha"); ylabel("CL"); title(["M=",Ms]); grid on
end