function [out1,out2,out3,out4] = Aerodynamics4Python(varP,npanels,angCal,Apexes,...
    varS,varA,M,supmob,deltaFs,deltaSs,dF,h,Neng,Qs,Nac,MTOM,hwinglet)%,...
% ...
% flaps,nflaps,slats,nslats,deltaFs,deltaSs,,elevator,deltaE)
%Aerodynamics esegue i calcoli aerodinamici completi per un velivolo usando
%le formule e i metodi del tool Aerodynamics_
%   varP: vettore che contiene su una riga le informazioni del pannello di
%       una superficie portante. [b,sweep,dihedral]
%   nPanels: vettore che contiene sulle colonne il numero di pannelli per
%       superficie portante,
%       [npanels(wing),npanel(horizontal),npanel(vertical)]'
%   Apexes: vettore contenente su una riga le coordinate X,Y,Z globali
%       dell'apice della rispettiva supeficie portante
%   varS: vettore contenente le informazioni geometriche dei profili
%   varA: vettore contenente le informazioni aerodinamiche dei profili a
%       dato Mach
%   angCal: vettore di angoli di calettamento
%   M: Mach di volo
%   supmob: matrice contenente le informazioni sulel usperfici mobili
%   deltaFs: deflessione dei flaps per date condizioni di volo
%   deltaSs: deflessione dei slats per date condizioni di volo
%   dF: diametro fusoliera
%   h: quota di volo
%   Neng: numero motori
%   Qs: fattori di interferenza, nell'ordine Qs(1) = QS (interferenza
%       singola superfice)
%       Qs(2) = QSN (interferenza gondola-superficie)
%   BarCoord: coordinate del baricentro del velivolo completo %AGGIUNGERE
%       NEL .py
%   Nac = [Ln, FF, Swet]: vettore contnente informazioni sulle nacelles
% Calcola il numero di superfici portanti dichiarate: nell'ordine devono
% essere Wing;Horizontal,Vertical,Canard
 
% % DEBUG
% PYTHON = 0;
% if PYTHON == 1
%     disp("save")
%     save("varDebug.mat","varP","npanels","angCal","Apexes",...
%         "varS","varA","M","supmob","deltaFs","deltaSs","dF","h",...
%         "Neng","Qs","Nac","MTOM","hwinglet")
% else
%     load("C:\Users\Ospitio\Documents\RCE\Tools\B747Tools\DEBUG\Copy_Aerodynamic_B747\matlab\varDebug.mat")
% end

%tic
nsurfs = 0; totpnls = 0;
for k=1:length(npanels)
    totpnls = totpnls + npanels(k);
    if npanels(k) ~= 0
        nsurfs = nsurfs + 1;
    end
end
% Calcola il numero di condizioni di volo
nFC = length(M);
flCnd = {'clean','Landing','Take-Off'};
Cd0Fus=NaN(nFC,1); CM0Fus = Cd0Fus; CmaFus = Cd0Fus; Cd0Nac = Cd0Fus; dCd0Und = Cd0Fus;
%% High Lift e Superfici Mobili
% Il vettore di input avrà sulla prima colonna l'ID della superficie a cui
% è associata la superficie mobile ( 1 ala, 2 piano orizz.), sulal seconda
% colonna avrà l'ID della superficie mobile (per l' ala 1=flap,2slat, per
% piano orizz. 1 = elevator). Le altre colonne saranno le grandezze delle
% superfici
nSm = length( supmob(:,1) ) - 1; app = NaN(1,nSm);
for l = 1:nSm
    %app è un vettore che contiene le differenze degli ID tra due posizioni
    %consecutive
    app(1,l) = supmob(l+1,1) - supmob(l,1) - 1;
end
app = [app,0]; app(app<1) = 0;%[app,zeros(2,1)]; app(app<1) = 0;
% Si vuoel creare una matrice che abbia tutte le superfici portanti, anche
% quelle su cui non sono definite le superfici mobili; nSm in questo caso
% serve a calcolare la dimensioen finale della matrice
nSm = nSm + 1 + sum( app(1,app(1,:)>0) );
c = NaN( nSm,length(supmob(1,:)) );
nSm = length( supmob(:,1) ) + 1;
k = 1; j = 1;
while j < nSm
    if app(1,j) < 1
        % Se in posizione j app(j) = 0 allora vuol dire che tra j e j+1
        % si sta ancora sulla stessa superficie
        c(k,:) = supmob(j,:);
        % d(k) = deltas(j);
        k = k +1; j = j + 1;
    else
        % Se invece la differenza tra gli ID è > 1 vuol dire che non è
        % stata assegnata almeno uan superficie intermedia
        c( k:k+app(1,j),: ) = [ (supmob(j,1):supmob(j,1)+app(1,j))', ...
            zeros( app(1,j)+1 ,length(supmob(1,:))-1 )] ;
        % d( k:k+app(1,j) ) = zeros(1,app(1,j)+1);
        k = k + app(1,j) + 1; j = j+1;
    end
end
% Se l'ID dell'iltima superficie è minore del numero di superfici
% allora riempie le righe con 0 ( due righe per superficie)
if c(end,1) < nsurfs
    nlast = length(c(:,1));
    c = [c; NaN( ( nsurfs - c(end,1) )*2 , length(supmob(1,:)) )];
    k = c(nlast,1); j = 1;
    while  k < nsurfs+1
        c(nlast+k:nlast+k+1,1) = ones(2,1)*( c(nlast,1)+j );
        k = k+2; j = j +1;
    end
end
nSm = length( c(:,1)  ); % due sezioni per superficie mobile
nmobsurf = zeros(nsurfs,1); % Contatore del numero di superfici mobili per superficie portante
hflag = repmat( {'-'},nSm,1);

k = 1; j = 1;
while k < nSm + 1
    switch c(k,1) % Scelta della superifice su cui è inserito
        case 1 % Ala
            switch c(k,2) % Identificazione del tipo di supmob
                case 1 %flap
                    hflag{j} = 'flaps';
                    nmobsurf( c(k,1) ) = nmobsurf( c(k,1) ) + 1;
                    k = k + 2;
                    % u = u + 1;
                case 2 %slat
                    hflag{j} = 'slats';
                    nmobsurf( c(k,1) ) = nmobsurf( c(k,1) ) + 1;
                    k = k + 2;
                    % uu = uu + 1;
                otherwise
                    k = k + 1;
            end
        case 2 % Piano Orizzontale
            switch c(k,2)
                case 1 % Elevator
                    hflag{j} = 'elevator';
                    nmobsurf( c(k,1) ) = nmobsurf( c(k,1) ) + 1;
                    k = k +2;
                otherwise
                    k = k + 1;
            end
        otherwise
            k = k + 1;
        % case 3 % Piano Verticale
        %     switch c(k,2)
        %     end
    end
            j = j+1;
end
nmobsurf(nmobsurf<1) = 1; % serve a fare in modo che nel ciclo si continui a scorrere hflag anche se la sup non ha sup. mob.
%% Calcoli Aerodinamici
% Preallocazione delle variabili di input
out1 = NaN(nsurfs,12);
out2 = NaN(totpnls,8);
out3 = NaN(nsurfs*nFC,16);
LiftingSurfs = WingClass.empty;

cnt = 0; cSt = 0; Neng = [Neng;zeros(nsurfs-1,1)];
%% Definizione delle Superfici e calcolo in condizioni pulite
for i = 1: nsurfs
    fuselage.dfmax = dF(i);
    % Definizione
    LiftingSurfs(i) = WingClass( varP(cnt+1:cnt+npanels(i),1),varP(cnt+1:cnt+npanels(i),2),...
        varP(cnt+1:cnt+npanels(i),3),angCal(i),Apexes(i,:),M,varS( cnt + i : cnt+npanels(i)+i,:), ...
        varA((cnt+i-1)*nFC + 1 : (cnt+npanels(i) + i)*nFC,:),...
        hflag( cSt+1 : cSt+nmobsurf(i) ),c( 2*cSt+1 : 2*(cSt+nmobsurf(i)),3:end ) );
    % Calcolo per ogni Mach
    LiftingSurfs(i).wing3Ddata(1) = LiftingSurfs(i).aero3Dwing(flCnd{1});

    % out2: dimensioni dei singoli pannlli delal superficie:
    % S, TR, xglob_root,yglob_root,zglob_root, xglob_tip,yglob_tip,zglob_tip
    for k=1:npanels(i)
        out2( cnt+k,: ) = [LiftingSurfs(i).panels(k).S,LiftingSurfs(i).panels(k).TR,...
            LiftingSurfs(i).panels(k).root.xglob,LiftingSurfs(i).panels(k).root.yglob,...
            LiftingSurfs(i).panels(k).root.zglob, LiftingSurfs(i).panels(k).tip.xglob,...
            LiftingSurfs(i).panels(k).tip.yglob,LiftingSurfs(i).panels(k).tip.zglob];
    end
    % Calcolo Drag in Condizioni Pulite
    m = 1;
    switch i
        case 1
            LiftingSurfs(1).hwglt = hwinglet;
        case 3
            % Nel piano verticale viene calcolata bv e Sv come se lil piano
            % fosse simmetrico (tipo un ala), quindi AR = (2bv)^2/(2Sv) = 2
            % AR_true. Non modifichiamo ancora S e b perchè
            % modificherebbero Sexp nel calcolo di Cd0
            LiftingSurfs(3).AR = LiftingSurfs(3).AR*0.5;
    end
    [LiftingSurfs(i).wing3Ddata(1).Cd0, LiftingSurfs(i).wing3Ddata(1).dCd0, ...
                LiftingSurfs(i).wing3Ddata.e] = LiftingSurfs(i).DragEstLifSur(...
                LiftingSurfs(i).wing3Ddata(1),fuselage,h,Neng(i),...
                Qs(nFC*(i-1)+m,1),Qs(nFC*(i-1)+m,2) );

    cnt = cnt + npanels(i); cSt = cSt + nmobsurf(i);
end
% Dopo aver calcolato il Cd0 del piano verticale si possono aggiornare Sv e
% bv -> Cd0 è proporzionale a Sexp/Sv
LiftingSurfs(3).Sw = LiftingSurfs(3).Sw*0.5; LiftingSurfs(3).bw = LiftingSurfs(3).bw * 0.5;

%% Calcoli di High-Lift per l'Ala
for m = 2:nFC
    LiftingSurfs(1).wing3Ddata(m) = LiftingSurfs(1).aero3Dwing(flCnd{m},M(m),deltaFs(m),deltaSs(m));
    % Calcolo dell'incremento di resistenza dovuto ai flaps
    [~, LiftingSurfs(1).wing3Ddata(m).dCd0, ...
        ~] = LiftingSurfs(1).DragEstLifSur(...
        LiftingSurfs(1).wing3Ddata(m),fuselage,h(m),Neng(1),...
        Qs(nFC*(1-1)+m,1),Qs(nFC*(1-1)+m,2),M(m));

end

%% Calcoli Downwash e Fusoliera
Kuc = [0,4.49e-5,3.16e-5];
for m = 1:nFC
    % Downwash
    [ LiftingSurfs(2).wing3Ddata(1).eps0(m),LiftingSurfs(2).wing3Ddata(1).depsda(m) ] = ...
        LiftingSurfs(2).downwashEstimation(...
        LiftingSurfs(1).meanprofile.xglob + LiftingSurfs(1).meanprofile.c*0.25 ,...
        LiftingSurfs(1).meanprofile.zglob,...
        LiftingSurfs(2).meanprofile.xglob + LiftingSurfs(2).meanprofile.c*0.25, ...
        LiftingSurfs(2).meanprofile.zglob,...
        LiftingSurfs(1).iAng, LiftingSurfs(1).wing3Ddata(m).alpha0l(m), ...
        LiftingSurfs(1).bw, LiftingSurfs(1).AR,LiftingSurfs(1).wing3Ddata(m).a(m),...
        LiftingSurfs(1).panels(1).sweep);
    % Fusoliera
    [dfmax,lf,~,hol075,dCdWindShield,db,Swet,CM0Fus(m),CmaFus(m)] = Fuselage( LiftingSurfs(1).wing3Ddata(m).alpha0l(m) - LiftingSurfs(1).iAng, LiftingSurfs(1).wing3Ddata(m).alpha0l(m), ...
        LiftingSurfs(1).wing3Ddata(m).alpha0l(m), LiftingSurfs(1).meanprofile.c,...
        LiftingSurfs(1).meanprofile.xglob + LiftingSurfs(1).meanprofile.c*0.25,...
        LiftingSurfs(1).panels(1).root.xglob, LiftingSurfs(1).panels(1).root.xglob + LiftingSurfs(1).panels(1).root.c, ...
        LiftingSurfs(1).wing3Ddata(m).a(m), LiftingSurfs(1).Sw,LiftingSurfs(1).panels(1).root.zglob, ...
        LiftingSurfs(2).wing3Ddata(1).depsda(m) );
    % Drag di Fusoliera
    [Cd0Fus(m), Cd0Nac(m)] = Drag_est_fus(M(m),h(m),lf,dfmax,4,db,...
        dCdWindShield,hol075,Swet,LiftingSurfs(1).Sw, Nac(1), Nac(2), Nac(3));
    % Undercarriage:
    dCd0Und(m) = MTOM^(1-0.215)*9.81*Kuc(m)/LiftingSurfs(1).Sw;
end

%% Scrittura degli Output

% Fusoliera + Nacelles + Undercarriages
out4 = [Cd0Fus(:), CM0Fus(:), CmaFus(:), Cd0Nac(:), dCd0Und(:),ones(nFC,1)*Swet ];
% Ala
i = 1;
for m = 1:nFC
    out3((i-1)*nFC+m,:) = [LiftingSurfs(i).wing3Ddata(m).a(m), LiftingSurfs(i).wing3Ddata(m).cl0(m),...
        LiftingSurfs(i).wing3Ddata(m).clstar(m),LiftingSurfs(i).wing3Ddata(m).clmax(m),...
        LiftingSurfs(i).wing3Ddata(m).alpha0l(m),LiftingSurfs(i).wing3Ddata(m).alphastar(m),...
        LiftingSurfs(i).wing3Ddata(m).alphamax(m),LiftingSurfs(i).wing3Ddata(m).cmac(m),...
        LiftingSurfs(i).wing3Ddata(m).Cd0(m),LiftingSurfs(i).wing3Ddata(m).dCd0(m),...
        LiftingSurfs(i).wing3Ddata(m).e(m),LiftingSurfs(i).wing3Ddata(m).Re(m),...
        LiftingSurfs(i).wing3Ddata(m).deltaF,LiftingSurfs(i).wing3Ddata(m).deltaS,0,0];
end
out1(i,:) = [LiftingSurfs(i).bw,LiftingSurfs(i).Sw,LiftingSurfs(i).TR,...
    LiftingSurfs(i).AR,LiftingSurfs(i).meanprofile.c,...
    LiftingSurfs(i).meanprofile.xglob,LiftingSurfs(i).meanprofile.yglob,...
    LiftingSurfs(i).meanprofile.zglob,LiftingSurfs(i).Swet,...
    LiftingSurfs(i).hwglt,LiftingSurfs(i).tau,LiftingSurfs(i).meanprofile.tc];
% Piano Orizzontale
i = 2;
out3((i-1)*nFC+1 : i*nFC,:) = [LiftingSurfs(i).wing3Ddata.a(:), LiftingSurfs(i).wing3Ddata.cl0(:),...
    LiftingSurfs(i).wing3Ddata.clstar(:),LiftingSurfs(i).wing3Ddata.clmax(:),...
    LiftingSurfs(i).wing3Ddata.alpha0l(:),LiftingSurfs(i).wing3Ddata.alphastar(:),...
    LiftingSurfs(i).wing3Ddata.alphamax(:),LiftingSurfs(i).wing3Ddata.cmac(:),...
    LiftingSurfs(i).wing3Ddata.Cd0(:),LiftingSurfs(i).wing3Ddata.dCd0(:),...
    LiftingSurfs(i).wing3Ddata.e(:),LiftingSurfs(i).wing3Ddata.Re(:),...
    zeros(nFC,4)];
out3((i-1)*nFC+1 : i*nFC ,end-1:end) = [ LiftingSurfs(i).wing3Ddata(1).eps0(:),LiftingSurfs(i).wing3Ddata(1).depsda(:)];
out1(i,:) = [LiftingSurfs(i).bw,LiftingSurfs(i).Sw,LiftingSurfs(i).TR,...
    LiftingSurfs(i).AR,LiftingSurfs(i).meanprofile.c,...
    LiftingSurfs(i).meanprofile.xglob,LiftingSurfs(i).meanprofile.yglob,...
    LiftingSurfs(i).meanprofile.zglob,LiftingSurfs(i).Swet,...
    LiftingSurfs(i).hwglt,LiftingSurfs(i).tau,LiftingSurfs(i).meanprofile.tc];
% Piano Verticale
i = 3;
out3((i-1)*nFC+1 : i*nFC,:) = [zeros(nFC,8),...
    LiftingSurfs(i).wing3Ddata(1).Cd0(:)*0.5 ,LiftingSurfs(i).wing3Ddata(1).dCd0(:),...
    zeros(nFC,1),LiftingSurfs(i).wing3Ddata(1).Re(:)*0.5,...
    zeros(nFC,4)];
out1(i,:) = [LiftingSurfs(i).bw,LiftingSurfs(i).Sw,LiftingSurfs(i).TR,...
    LiftingSurfs(i).AR,LiftingSurfs(i).meanprofile.c,...
    LiftingSurfs(i).meanprofile.xglob,LiftingSurfs(i).meanprofile.yglob,...
    LiftingSurfs(i).meanprofile.zglob,LiftingSurfs(i).Swet,...
    LiftingSurfs(i).hwglt,LiftingSurfs(i).tau,LiftingSurfs(i).meanprofile.tc];





for k=1:npanels(i)
    out2( cnt+k,: ) = [LiftingSurfs(i).panels(k).S,LiftingSurfs(i).panels(k).TR,...
        LiftingSurfs(i).panels(k).root.xglob,LiftingSurfs(i).panels(k).root.yglob,...
        LiftingSurfs(i).panels(k).root.zglob, LiftingSurfs(i).panels(k).tip.xglob,...
        LiftingSurfs(i).panels(k).tip.yglob,LiftingSurfs(i).panels(k).tip.zglob];
end





% %     % b, S, TR, AR, c_mac, xglob_mac, yglob_mac, zglob_mac, Swet, hwglt,
% %     %   tau, tc_averange
% %     out1(i,:) = [LiftingSurfs(i).bw,LiftingSurfs(i).Sw,LiftingSurfs(i).TR,...
% %         LiftingSurfs(i).AR,LiftingSurfs(i).meanprofile.c,...
% %         LiftingSurfs(i).meanprofile.xglob,LiftingSurfs(i).meanprofile.yglob,...
% %         LiftingSurfs(i).meanprofile.zglob,LiftingSurfs(i).Swet,...
% %         LiftingSurfs(i).hwglt,LiftingSurfs(i).tau,LiftingSurfs(i).meanprofile.tc];





% % 
% % %% 
% % for i = 1:nsurfs
% %     % Definizione del Diametro di Fusoliera per l'iesimo componente
% %     fuselage.dfmax = dF(i);
% %     % Definizione della SUperficie i-esima (ala, piano di coda ...)
% %     LiftingSurfs(i) = WingClass( varP(cnt+1:cnt+npanels(i),1),varP(cnt+1:cnt+npanels(i),2),...
% %         varP(cnt+1:cnt+npanels(i),3),angCal(i),Apexes(i,:),M,varS( cnt + i : cnt+npanels(i)+i,:), ...
% %         varA((cnt+i-1)*nFC + 1 : (cnt+npanels(i) + i)*nFC,:),...
% %         hflag( cSt+1 : cSt+nmobsurf(i) ),c( 2*cSt+1 : 2*(cSt+nmobsurf(i)),3:end ) );
% %     % Calcola le grandezze dell'ala pulita: sono escluse dal ciclo su m
% %     % perchè in questo modo calcola in una riga le grandezze 3D per tutte
% %     % le condizioni
% %     % m = 1;
% %     LiftingSurfs(i).wing3Ddata(1) = LiftingSurfs(i).aero3Dwing(flCnd{1});
% %     if i == 1
% %         % CALCOLI ALA
% %         LiftingSurfs(i).hwglt = hwinglet;
% % 
% %         % Con l'ala è necessario cal
% %         for m = 1:nFC
% %             LiftingSurfs(i).wing3Ddata(m) = LiftingSurfs(i).aero3Dwing(flCnd{m},M(m),deltaFs(m),deltaSs(m));
% % 
% %             [LiftingSurfs(i).wing3Ddata(m).Cd0, LiftingSurfs(i).wing3Ddata(m).dCd0, ...
% %                 LiftingSurfs(i).wing3Ddata(m).e] = LiftingSurfs(i).DragEstLifSur(...
% %                 LiftingSurfs(i).wing3Ddata(m),fuselage,h(m),Neng,...
% %                 Qs(nFC*(i-1)+m,1),Qs(nFC*(i-1)+m,2),M(m));
% %             out3((i-1)*nFC+m,:) = [LiftingSurfs(i).wing3Ddata(m).a(m), LiftingSurfs(i).wing3Ddata(m).cl0(m),...
% %                 LiftingSurfs(i).wing3Ddata(m).clstar(m),LiftingSurfs(i).wing3Ddata(m).clmax(m),...
% %                 LiftingSurfs(i).wing3Ddata(m).alpha0l(m),LiftingSurfs(i).wing3Ddata(m).alphastar(m),...
% %                 LiftingSurfs(i).wing3Ddata(m).alphamax(m),LiftingSurfs(i).wing3Ddata(m).cmac(m),...
% %                 LiftingSurfs(i).wing3Ddata(m).Cd0(m),LiftingSurfs(i).wing3Ddata(m).dCd0(m),...
% %                 LiftingSurfs(i).wing3Ddata(m).e(m),LiftingSurfs(i).wing3Ddata(m).Re(m),...
% %                 LiftingSurfs(i).wing3Ddata(m).deltaF,LiftingSurfs(i).wing3Ddata(m).deltaS,0,0];
% % 
% %         end
% %     else
% %         % Neng = 0 se la superficie non è l'ala
% %         m = 1;
% %         %LiftingSurfs(i).wing3Ddata(m) = LiftingSurfs(i).aero3Dwing(flCnd{1});
% %         if i == 3
% %             % PIANO VERTICALE: servono solo i dati di resistenza che vanno
% %             % divisi per 2 poichè viene calcolata una Swet che è il doppio
% %             % di quella effettiva
% %             LiftingSurfs(i).AR = LiftingSurfs(i).AR*0.5; LiftingSurfs(i).bw = LiftingSurfs(i).bw*0.5;
% %             LiftingSurfs(i).Sw = LiftingSurfs(i).Sw*0.5;
% %             % Calcolo della resistenza: ATTENZIONE dal momento che usa la
% %             % funzione dell' oggetto superficie i-esima Sw sarà quello della
% %             % superficie in questione e non dell'ala, quindi il Cd0 è riferito
% %             % alla superficie i-esima e non a Sw
% %             [LiftingSurfs(i).wing3Ddata.Cd0, LiftingSurfs(i).wing3Ddata.dCd0, ...
% %                 LiftingSurfs(i).wing3Ddata.e] = LiftingSurfs(i).DragEstLifSur(...
% %                 LiftingSurfs(i).wing3Ddata,fuselage,h,0,...
% %                 Qs(nFC*(i-1)+m,1),Qs(nFC*(i-1)+m,2) );
% %             out3((i-1)*nFC+1 : i*nFC,:) = [zeros(nFC,8),...
% %                 LiftingSurfs(i).wing3Ddata(1).Cd0(:)*0.5 ,LiftingSurfs(i).wing3Ddata(1).dCd0(:),...
% %                 zeros(nFC,1),LiftingSurfs(i).wing3Ddata(1).Re(:)*0.5,...
% %                 zeros(nFC,4)];
% % 
% %         else
% %             % Calcolo della resistenza: ATTENZIONE dal momento che usa la
% %             % funzione dell' oggetto superficie i-esima Sw sarà quello della
% %             % superficie in questione e non dell'ala, quindi il Cd0 è riferito
% %             % alla superficie i-esima e non a Sw
% %             [LiftingSurfs(i).wing3Ddata.Cd0, LiftingSurfs(i).wing3Ddata.dCd0, ...
% %                 LiftingSurfs(i).wing3Ddata.e] = LiftingSurfs(i).DragEstLifSur(...
% %                 LiftingSurfs(i).wing3Ddata,fuselage,h,0,...
% %                 Qs(nFC*(i-1)+m,1),Qs(nFC*(i-1)+m,2) );
% %             out3((i-1)*nFC+1 : i*nFC,:) = [LiftingSurfs(i).wing3Ddata(1).a(:), LiftingSurfs(i).wing3Ddata(1).cl0(:),...
% %                 LiftingSurfs(i).wing3Ddata(1).clstar(:),LiftingSurfs(i).wing3Ddata(1).clmax(:),...
% %                 LiftingSurfs(i).wing3Ddata(1).alpha0l(:),LiftingSurfs(i).wing3Ddata(1).alphastar(:),...
% %                 LiftingSurfs(i).wing3Ddata(1).alphamax(:),LiftingSurfs(i).wing3Ddata(1).cmac(:),...
% %                 LiftingSurfs(i).wing3Ddata(1).Cd0(:),LiftingSurfs(i).wing3Ddata(1).dCd0(:),...
% %                 LiftingSurfs(i).wing3Ddata(1).e(:),LiftingSurfs(i).wing3Ddata(1).Re(:),...
% %                 zeros(nFC,4)];
% %         end
% %         if i == 2
% %             % PIANO ORIZZONTALE: calcolo downwash
% %             for m = 1:nFC
% %                 [ LiftingSurfs(i).wing3Ddata(1).eps0(m),LiftingSurfs(i).wing3Ddata(1).depsda(m) ] = ...
% %                     LiftingSurfs(i).downwashEstimation(...
% %                     LiftingSurfs(1).meanprofile.xglob + LiftingSurfs(1).meanprofile.c*0.25 ,...
% %                     LiftingSurfs(1).meanprofile.zglob,...
% %                     LiftingSurfs(i).meanprofile.xglob + LiftingSurfs(i).meanprofile.c*0.25, ...
% %                     LiftingSurfs(i).meanprofile.zglob,...
% %                     LiftingSurfs(i-1).iAng, LiftingSurfs(i-1).wing3Ddata(m).alpha0l(m), ...
% %                     LiftingSurfs(i-1).bw, LiftingSurfs(i-1).AR,LiftingSurfs(i-1).wing3Ddata(m).a(m),...
% %                     LiftingSurfs(i-1).panels(1).sweep);
% %             end
% %         else
% %             % Downwash nullo
% %             LiftingSurfs(i).wing3Ddata(1).eps0 = zeros(nFC,1); LiftingSurfs(i).wing3Ddata(1).depsda = zeros(nFC,1);
% %         end
% %         out3((i-1)*nFC+1 : i*nFC ,end-1:end) = [ LiftingSurfs(i).wing3Ddata(1).eps0(:),LiftingSurfs(i).wing3Ddata(1).depsda(:)];
% %     end
% %     %end
% %     %end
% %     % Scrittura degli output
% %     % b, S, TR, AR, c_mac, xglob_mac, yglob_mac, zglob_mac, Swet, hwglt,
% %     %   tau, tc_averange
% %     out1(i,:) = [LiftingSurfs(i).bw,LiftingSurfs(i).Sw,LiftingSurfs(i).TR,...
% %         LiftingSurfs(i).AR,LiftingSurfs(i).meanprofile.c,...
% %         LiftingSurfs(i).meanprofile.xglob,LiftingSurfs(i).meanprofile.yglob,...
% %         LiftingSurfs(i).meanprofile.zglob,LiftingSurfs(i).Swet,...
% %         LiftingSurfs(i).hwglt,LiftingSurfs(i).tau,LiftingSurfs(i).meanprofile.tc];
% %     % out2: dimensioni dei singoli pannlli delal superficie:
% %     % S, TR, xglob_root,yglob_root,zglob_root, xglob_tip,yglob_tip,zglob_tip
% %     for k=1:npanels(i)
% %         out2( cnt+k,: ) = [LiftingSurfs(i).panels(k).S,LiftingSurfs(i).panels(k).TR,...
% %             LiftingSurfs(i).panels(k).root.xglob,LiftingSurfs(i).panels(k).root.yglob,...
% %             LiftingSurfs(i).panels(k).root.zglob, LiftingSurfs(i).panels(k).tip.xglob,...
% %             LiftingSurfs(i).panels(k).tip.yglob,LiftingSurfs(i).panels(k).tip.zglob];
% %     end
% %     cnt = cnt + npanels(i); cSt = cSt + nmobsurf(i);
% % end
% % Kuc = [0,4.49e-5,3.16e-5];
% % for m = 1:nFC
% %     % Calcolo Grandezze Fusoliera:
% %     % dfmax,lf,Xcabin_start,hol075,Delta_CD_wshield,db,Swet,CM0_Fuso,CMa_Fuso,CDups_Fuso
% %     [dfmax,lf,~,hol075,dCdWindShield,db,Swet,CM0Fus(m),CmaFus(m)] = Fuselage( LiftingSurfs(1).wing3Ddata(m).alpha0l(m) - LiftingSurfs(1).iAng, LiftingSurfs(1).wing3Ddata(m).alpha0l(m), ...
% %         LiftingSurfs(1).wing3Ddata(m).alpha0l(m), LiftingSurfs(1).meanprofile.c,...
% %         LiftingSurfs(1).meanprofile.xglob + LiftingSurfs(1).meanprofile.c*0.25,...
% %         LiftingSurfs(1).panels(1).root.xglob, LiftingSurfs(1).panels(1).root.xglob + LiftingSurfs(1).panels(1).root.c, ...
% %         LiftingSurfs(1).wing3Ddata(m).a(m), LiftingSurfs(1).Sw,LiftingSurfs(1).panels(1).root.zglob, ...
% %         LiftingSurfs(2).wing3Ddata(1).depsda(m) );
% %     % M,h,Lf,df,Type_of_surface_f,db,Dcd_windshield,h_l_075,Swet
% %     [Cd0Fus(m), Cd0Nac(m)] = Drag_est_fus(M(m),h(m),lf,dfmax,4,db,...
% %         dCdWindShield,hol075,Swet,LiftingSurfs(1).Sw, Nac(1), Nac(2), Nac(3));
% %     % Undercarriage:
% %     dCd0Und(m) = MTOM^(1-0.215)*9.81*Kuc(m)/LiftingSurfs(1).Sw;
% % end
% % % PIANO VERTICALE: dimezza b,S,AR
% % %out1(3,1) = out1(3,1) *0.5; out1(3,2) = out1(3,2) *0.5; out1(3,3) = out1(3,3) *0.5;
% % 
% % out4 = [Cd0Fus(:), CM0Fus(:), CmaFus(:), Cd0Nac(:), dCd0Und(:),ones(nFC,1)*Swet ];
% % %toc
% % end