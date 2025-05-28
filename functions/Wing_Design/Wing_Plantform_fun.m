function [f,ax,yroot,ykink,ytip,croot,ckink,ctip,xLE_root,xLE_kink, ...
    xLE_tip,A1_Sw] = Wing_Plantform_fun( Sw,bw,sweep,tr1,...
    ykink_b,tr2,A1_Sw )
%UNTITLED Summary of this function goes here
%INPUT
% Sw        : superficie alare
% ARw       : aspect ratio dell'ala
% sweep     : angolo di freccia Leading Edge
% tr1       : taper ratio c_kink/c_root (panel 1)
% ykink_b   : posizione del kink lungo lo span
% tr2       : taper ratio c_tip/c_root (panel 2)
% eps_root  : twist del profilo di root
% eps_tip   : twist del profilo di tip
% A1_Sw     : percentuatel dl'area del panel 1 rispetto alla totale
%OUTPUT
% MAC       : corda media aerodinamica
% k1        : coeffiente per il triangolo 1 dell'ala
% k2        : coeffiente per il triangolo 2 dell'ala
% k3        : coeffiente per il triangolo 3 dell'ala
% yroot     : posizione della root lungo lo span (0)
% ykink     : posizione della kink lungo lo span
% ytip      : posizione della tip  lungo lo span
% croot     : corda del profilo di root
% ckink     : corda del profilo di kink
% ctip      : corda del profilo di tip
% xLE_root  : posizione del LE della root
% xLE_kink  : posizione del LE della kink
% xLE_tip   : posizione del LE della tip
yroot    = 0;    
ykink    = ykink_b*bw/2;      
ytip     = bw/2;
bkink    = bw*ykink_b/2; 
tol = 1e-2; diff = 1;
%% Graphics
f = figure();
f.Name = 'Wing Plantform Parameters';
% adding plots
ax = axes('Parent',f); ax.Units = 'pixels';
ax.Position = [75 75 325 280]; cnt = 1;
%% Wing Area Convercenge Loop
while ( diff > tol ) && A1_Sw > 0
    A1d      = A1_Sw*Sw;
    A        = [bkink bkink;tr1 -1];
    noti     = [2*A1d ; 0];
    Ainv     = inv(A);
    croot    = Ainv(1,:)*noti;
    ckink    = Ainv(2,:)*noti;
    ctip     = tr2*croot;
    xLE_root = 0;                                xTE_root = xLE_root + croot;
    xLE_kink = xLE_root+tand(sweep)*ykink;       xTE_kink = xLE_kink + ckink;
    xLE_tip  = xLE_root+tand(sweep)*ytip;        xTE_tip  = xLE_tip + ctip;
    % Wing Plot
    wingPlot(f,ax,croot,ckink,ctip,xLE_root,xLE_kink,xLE_tip,yroot,ykink,ytip);
    % Wing Aerea Check
    A1   = .5*(ckink+croot)*(ykink-yroot);
    A2   = .5*(ckink+ctip)*(ytip-ykink);
    Atot = 2*(A1+A2);   diff = 1+abs( 1-Atot/Sw );
    fprintf("La differenza tra Sw scelta nel Sizing e la Sw calcolata con il" + ...
        "metodo delle aeree è: %.3f\n",diff*100);
    if cnt > 1 || diff > tol
        disp('Change What')
        disp(['1. A1/Sw: Current ',num2str(A1_Sw)])
        disp(['2. TRi  : Current ',num2str(tr2)])
        chg = input('>>');
        switch chg
            case 1
                disp('Insert new Ai/Sw')
                A1_Sw = input('>>');
            case 2
                disp('Insert new TR')
                tr2 = input('>>');
            otherwise
                break
        end
        
    end
    cnt = cnt + 1;
end

end

function wingPlot(f,ax,croot,ckink,ctip,xLE_root,xLE_kink,xLE_tip,yroot,ykink,ytip)
    pl = plot(ax,[xLE_root,xLE_kink,xLE_tip],[yroot,ykink,ytip]); hold on
    p2 = plot(ax,[xLE_root,xLE_kink,xLE_tip]+[croot,ckink,ctip],[yroot,ykink,ytip]);
    p3 = plot(ax,[xLE_tip,xLE_tip]+[0,ctip],[ytip,ytip]); hold off
end