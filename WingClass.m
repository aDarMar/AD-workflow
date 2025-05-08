classdef WingClass %< handle %<WingClass è una sottoclassed della classe predefinita handle
    properties
        % Geometry
            % Wing
        bw {mustBeNumeric, mustBePositive}              %apertura [m]
        Sw {mustBeNumeric, mustBePositive}              %Superficie (m^2)
        TR {mustBeNumeric, mustBePositive}              %Taper Ratio [-]
        AR {mustBeNumeric, mustBePositive}
        Swet {mustBeNumeric, mustBePositive}            %Area bagnata [m^2]
        iAng {mustBeNumeric}                             %Angolo di Calettamento dells Superficie
        hwglt {mustBeNumeric}                           %altezza winglet [m]
        npanels
        panels
            % High-Lifting
        nflaps {mustBeNumeric} 
        flaps
        cfocAvg {mustBeNumeric} %cf/c medio
        Sflaps {mustBeNumeric, mustBePositive} %Area Totale dei Flaps
        nslats {mustBeNumeric}
        Sslats {mustBeNumeric, mustBePositive} %Area Totale degli Slats
        slats
        csocAvg {mustBeNumeric}  %cs/c medio
            % Elevator
            commandSurf
        % Dati Intermedi
        meanprofile
        tau {mustBeNumeric, mustBePositive} % rAPPORTO SPESSORI
        % Dati per Calcolo Resistenza
        
        liftSurfFlag % flag che specifica il tipo di superficie
        % Aerodynamic Data
        wing3Ddata
        prf3DClean
        % Dati Aggiuntivi Piano di Coda 
        eps0 {mustBeNumeric}    %angolo di downwash per alpha = 0
        depsda {mustBeNumeric}  %Gradiente di downwash
        % High Lift Devices
        prf3Dflaps
        prf3Dflapslat

        % Global Coordinates 
        Xapex {mustBeNumeric}
        Yapex {mustBeNumeric}
        Zapex {mustBeNumeric}
        TST
    end
    methods
        function obj = WingClass(bs,sweeps,dihedrals,iang,apexC,M,...
                sectsGeom, sectsAero,...
                HLflag,sectsHL,...
                misc) %Costruttore
            if nargin == 1
                % Lettura XML
            else
                obj.TST = 1;
                nM = length(M); % Numero di condizioni di volo
                % Assegnazione Coordinate Apice Superficie
                obj.Xapex = apexC(1);
                obj.Yapex = apexC(2);
                obj.Zapex = apexC(3);
                % Assegnazione Angolo di Calettamento
                obj.iAng = iang;
                % Inizializzazione Pannelli
                obj.meanprofile = ProfileClass.empty;
                obj.npanels     = length(bs);
                obj.panels      = PanelClass.empty; % Costruisce un array di oggetti panels
                for i=1:obj.npanels
                    obj.panels(i) = PanelClass(bs(i),sweeps(i),dihedrals(i),M,...
                        sectsGeom(i ,:),sectsAero((i-1)*nM+1 : i*nM,:),...
                        sectsGeom(i+1,:),sectsAero(i*nM+1 : (i+1)*nM,:));
                end
                % Inizializzazione dati Winglet e Suèerficie Bagnata
                obj.Swet = 1; obj.hwglt = 0;
                if nargin == 11
                    obj.hwglt = misc(1);
                end
            end
            % Inizializzazione dei Dati di Downwash
            obj.eps0 = NaN; obj.depsda = NaN; obj.tau = 1;
            % Assegnazione delle coordinate ai profili
            i = 1;
            obj.panels(i).globalCoords( obj.Xapex,obj.Yapex,obj.Zapex ); 
            obj.panels(i).
            for i=2:obj.npanels
                obj.panels(i).globalCoords( obj.panels(i-1).tip.xglob,...
                    obj.panels(i-1).tip.yglob,obj.panels(i-1).tip.zglob )                
            end

            % Calcolo Grandezze Ala
            [obj.bw, obj.Sw] = geomCalc(obj);
            obj.TR = obj.panels(end).tip.c/obj.panels(1).root.c;
            obj.AR = obj.bw^2/obj.Sw;
            
            % High Lift ed Equilibratore
            obj.nflaps = 0; obj.nslats = 0;
            obj.cfocAvg = NaN; obj.csocAvg = NaN; 
            obj.Sflaps = 0.01; obj.Sslats = 0.01;
            % Definizione di Slat e Flap: sono degli oggetti della classe
            % Panel
            if exist('HLflag',"var")
                obj.flaps = PanelClass.empty;
                obj.slats = PanelClass.empty;
                obj.commandSurf = PanelClass.empty;
                sectsHL = [sectsHL(:,1)*0.5*obj.bw,sectsHL]; % Trasforma le coordinate adimensionali in dimensionali
                m = 1; mm = 1;
                while mm < length(HLflag) + 1
                    switch HLflag{mm}
                        case 'flaps'
                            obj.nflaps = obj.nflaps + 1;
                            % Dati Geometrici delle sezioni ei profili
                            varG = NaN(2,8+3);
                            % 2 coordinate di flap e 8 valori
                            varG(:,9:11) = sectsHL(m:m+1,2:end); %copia le grandezze di input dei flaps nella var temp
                            varG(:,9) = varG(:,9)*0.5*obj.bw;   %trasforma dy/b in coord. dimensionali
                            % Dati Aerodinamici delle Sezioni dei Profili
                            [varG,varA,coords] = HLAssign(obj,varG,nM);
                            obj.flaps(obj.nflaps) = PanelClass(varG(2,9) - varG(1,9),0,0,M,...
                                varG(1,1:8),varA(1:nM,:),varG(2,1:8),varA(nM+1:2*nM,:),HLflag{mm},sectsHL(m:m+1,:));
                            % Assegna le coordinate globali
                            obj.flaps(obj.nflaps).root.xglob = coords(1,1); obj.flaps(obj.nflaps).root.zglob = coords(1,2);
                            obj.flaps(obj.nflaps).tip.xglob = coords(2,1); obj.flaps(obj.nflaps).tip.zglob = coords(2,2);
                            m = m + 2; mm = mm + 1;
                        case 'slats'
                            obj.nslats = obj.nslats + 1;
                            % Dati Geometrici delle sezioni ei profili
                            varG = NaN(2,8+3);
                            % 2 coordinate di flap e 8 valori
                            varG(:,9:11) = sectsHL(m:m+1,2:end); %copia le grandezze di input dei flaps nella var temp
                            varG(:,9) = varG(:,9)*0.5*obj.bw;   %trasforma dy/b in coord. dimensionali
                            % Dati Aerodinamici delle Sezioni dei Profili
                            [varG,varA,coords] = HLAssign(obj,varG,nM);
                            obj.slats(obj.nslats) = PanelClass(varG(2,9) - varG(1,9),0,0,M,...
                                varG(1,1:8),varA(1:nM,:),varG(2,1:8),varA(nM+1:2*nM,:),HLflag{mm},sectsHL(m:m+1,:));
                            % Assegna le coordinate globali
                            obj.slats(obj.nslats).root.xglob = coords(1,1); obj.slats(obj.nslats).root.zglob = coords(1,2);
                            obj.slats(obj.nslats).tip.xglob = coords(2,1); obj.slats(obj.nslats).tip.zglob = coords(2,2);
                            m = m + 2; mm = mm + 1;
                        case 'elevator'
                            % Provvisorio, da finire quando si avranno le
                            % informazioni sull'elevator
                            obj.commandSurf = PanelClass(1,0,0,M,...
                                ones(1,8),ones(nM,8),ones(1,8),ones(nM,8),...
                                HLflag{mm},sectsHL(m:m+1,:));
                            m = m + 2; mm = mm + 1;
                        otherwise
                            mm = mm + 1;
                    end
                end
            end

            %Creazione del Profilo Medio
            [tmp1,tmp2] = meanProfileMod(obj);
            obj.meanprofile = ProfileClass(tmp1,tmp2,M);
            [~,xm,ym,zm] = macCalc(obj);
            obj.meanprofile.xglob = xm;
            obj.meanprofile.yglob = ym;
            obj.meanprofile.zglob = zm;


        end

        function Sw2 = sweepChange(~,Sw1,c1,c2,AR,TR)
            %sweepChange cambia l'angolo di freccia
            %   Sw1, c1 angolo di freccia (in gradi) alla percentuale c1 di corda. c2 percentuale
            %   di corda a cui si vuole calcolare lo sweep. A superficie dell'ala, TR
            %   taper ratio
            Sw1 = Sw1*pi/180; %Porta in radianti
            if (c1>1 || c1<0) || (c2>1 || c2<0)
                error('Out of bounds')
            end
            Sw2 = tan(Sw1) - 4/AR*(c2 - c1)*(1 - TR)/(1+TR);
            Sw2 = atan(Sw2);
            Sw2 = Sw2*180/pi; %[deg]
        end

        function [b,S] = geomCalc(obj)
            b = 0; S = 0;
            for i=1:obj.npanels
                b = b + obj.panels(i).b;
                S = S + obj.panels(i).S;
            end
            b = b*2; S = S*2;
        end

        function [mac, xmac, ymac, zmac] = macCalc(obj)
            mac = 0;
            xmac = 0;
            ymac = 0;
            zmac = 0;
            for i = 1:obj.npanels
                mac = mac + obj.panels(i).mac * 2*obj.panels(i).S/obj.Sw; %il 2 sta perchè panel.S è l'area del pannello e wing.S è l'area totale dell'ala
                xmac = xmac + (obj.panels(i).xmac + obj.panels(i).root.xglob) * 2*obj.panels(i).S/obj.Sw;
                % CONTROLLARE: sono ottenuti come medie pesate, vedere se è
                % vero.
                ymac = ymac + (obj.panels(i).ymac + obj.panels(i).root.yglob) * 2*obj.panels(i).S/obj.Sw; % CONTROLLARE!
                zmac = zmac + (obj.panels(i).zmac + obj.panels(i).root.zglob) * 2*obj.panels(i).S/obj.Sw;
            end
        end

        function avg = weightAvg(obj,grand)
            % Funzione che calcola i Ki delle medie pesate con le corde
            i = 1;
            avg = grand(:,i).*2.*obj.panels(i).b.*obj.panels(i).root.c.*0.5./obj.Sw;
            for i=1:obj.npanels-1
                avg = avg + grand(:,i+1).*2*(obj.panels(i).b + obj.panels(i+1).b)*obj.panels(i+1).root.c*0.5/obj.Sw;
                %                K(i) = (obj.panels(i-1).b + obj.panels(i).b)*obj.panels(i).root.c*0.5;
                %                avg = avg * grand(i)*2*K(i)/obj.Sw;
            end
            i = obj.npanels;
            avg = avg + grand(:,i+1).*2*obj.panels(i).b.*obj.panels(i).tip.c*0.5./obj.Sw;
        end
        
        function [vout,vout2] = meanProfileMod(obj)
            %meanProfile: calcola le grandezze medie dell'ala.
            % Funzione che calcola il profilo medio
            av = NaN(1,obj.npanels+1);
            dYv = av;xrtcv = av; xtrUpv = av; xrtLowv = av; tcv = av;
            nM = length(obj.panels(1).root.a);
            av = NaN(nM,obj.npanels+1);
            alphamaxv = av; alpha0lv = av; cmacv = av; alphastarv = av;
            cl0v = av; clstarv = av; clmaxv = av;
            
            for i = 1:obj.npanels
                tcv(i) = obj.panels(i).root.tc;
                xrtcv(i) = obj.panels(i).root.xtc;
                xtrUpv(i) = obj.panels(i).root.xtrUp;
                xrtLowv(i) = obj.panels(i).root.xtrLow;
                dYv(i) = obj.panels(i).root.dY;
                
                av(:,i) = obj.panels(i).root.a(:);
                cl0v(:,i) = obj.panels(i).root.cl0(:);
                clstarv(:,i) = obj.panels(i).root.clstar(:);
                clmaxv(:,i) = obj.panels(i).root.clmax(:);
                alphamaxv(:,i) = obj.panels(i).root.alphamax - obj.panels(i).root.eps(:);
                alpha0lv(:,i) = obj.panels(i).root.alpha0l - obj.panels(i).root.eps(:);
                alphastarv(:,i) = (obj.panels(i).root.clstar(:) -  obj.panels(i).root.cl0(:))./obj.panels(i).root.a(:);
                cmacv(:,i) = obj.panels(i).root.cmac(:);
            end
            i = obj.npanels;
            tcv(i+1) = obj.panels(i).tip.tc;
            xrtcv(i+1) = obj.panels(i).tip.xtc;
            xtrUpv(i+1) = obj.panels(i).tip.xtrUp;
            xrtLowv(i+1) = obj.panels(i).tip.xtrLow;
            dYv(i+1) = obj.panels(i).tip.dY;

            av(:,i+1) = obj.panels(i).tip.a;
            cl0v(:,i+1) = obj.panels(i).tip.cl0;
            clstarv(:,i+1) = obj.panels(i).tip.clstar;
            clmaxv(:,i+1) = obj.panels(i).tip.clmax;
            alphamaxv(:,i+1) = obj.panels(i).tip.alphamax - obj.panels(i).tip.eps;
            alpha0lv(:,i+1) = obj.panels(i).tip.alpha0l - obj.panels(i).tip.eps;
            alphastarv(:,i+1) = (obj.panels(i).tip.clstar -  obj.panels(i).tip.cl0)./obj.panels(i).tip.a;
            cmacv(:,i+1) = obj.panels(i).tip.cmac;

            [cm,~,~,~] = macCalc(obj);
            vout = [cm,1,obj.weightAvg(tcv),1,obj.weightAvg(xrtcv),obj.weightAvg(xtrUpv),...
                obj.weightAvg(xrtLowv),obj.weightAvg(dYv)];
            vout2 = [obj.weightAvg(av),...
                obj.weightAvg(cl0v),obj.weightAvg(clstarv),obj.weightAvg(clmaxv),...
                obj.weightAvg(alphamaxv),obj.weightAvg(alpha0lv),...
                obj.weightAvg(alphastarv),obj.weightAvg(cmacv)];

        end
        
        function [varG,varA,coords] = HLAssign(obj,varG,nM)
            %HLAssign: ricava per interpolazione i valori delle
            %caratteristiche geometriche e aerodinamiche dei profili che
            %delimitano i flaps, cui posizione sull'apertura è contenuta
            %nella colonna 9 di varG
            %   varG: vettore contenente le dimensioni geometriche degli
            %       ipersostentatori
            %   nM: numero di condizioni di volo
            %   varA: vettore contenente le caratteristiche aerodinamiche
            %       dei profili che delimitano i flaps
            varA = NaN(2*nM,8);
            %nM = length(varA(:,1))*0.5;
            ny = 0; j=1;
            coords = NaN(2,2);
            while ny<2
                [varGtmp,ypts] = obj.panels(j).panelInterp(...
                    varG(:,9));

                if ~isempty(ypts) && ny == 0
                    varA(1:length(ypts)*nM,:) = obj.panels(j).panelInterp(varG(:,9),8); % a
                    varG(1:length(ypts),1:8) = varGtmp;
                    coords(1:length(ypts),1) = obj.panels(j).root.xglob +  (varG(1:length(ypts),9) - obj.panels(j).root.yglob)*tan( obj.panels(j).sweep*pi/180); %coordinata x
                    coords(1:length(ypts),2) = obj.panels(j).root.zglob +  (varG(1:length(ypts),9) - obj.panels(j).root.yglob)*tan( obj.panels(j).dihedral*pi/180); %coordinata z
                elseif ~isempty(ypts) %&& ny == 1
                    varA(nM+1:2*nM,:) = obj.panels(j).panelInterp(varG(:,9),8); % a
                    varG(2,1:8) = varGtmp;
                    coords(2,1) = obj.panels(j).root.xglob +  (varG(2,9) - obj.panels(j).root.yglob)*tan( obj.panels(j).sweep*pi/180); %coordinata x
                    coords(2,2) = obj.panels(j).root.zglob +  (varG(2,9) - obj.panels(j).root.yglob)*tan( obj.panels(j).dihedral*pi/180); %coordinata z
                end
                ny = ny+length(ypts);
                j = j+1;
            end

        end
        
        function obc = tridProfInit(obj)
            %tridProfInit: funzione che inizializza un oggetto della classe
            %profilo per renderlo adattop al calcolo delle grandezze 2D
            %dell'ala da parte di aero3Dwing
            %   0bc: oggetto classe profile che deve essere inizializzato
            [tmp1,tmp2] = meanProfileMod(obj);
            
            % if nargin == 1
            %M = obj.meanprofile.M;
            % end
            obc = ProfileClass(tmp1,tmp2,obj.meanprofile.M);
            [~,xm,ym,zm] = macCalc(obj);
            obc.xglob = xm;
            obc.yglob = ym;
            obc.zglob = zm;
        end
        %% Portanza 3D

        function profClass = aero3Dwing(obj, flg, M, deltaF, deltaS,Fcalc)
            %aero3Dwing calcola le caratteristiche aerodinamiche dell'ala
            %tridimensionale
            %   profClass oggetto classe profilo sul quale salvare i dati
            %       aerodinamici
            %   Fcalc: se definito forza il calcolo dei coefficienti
            %   dell'ala pulita nel calcolo dell'ala con ipersost 
            if nargin > 2
                % Controlla se viene assegnato il Mach, altrimenti calcola
                % i dati 3D per ogni mach immagazzinato nei profili
                
                % Controlla se nell'oggetto sono salvati i dati
                % aerodinamici al amch imposto e ritorna un vettore di
                % indici delle posizioni a cui corrispondono i dati
                % aerodinamici al mach fissato
                nM = length(obj.meanprofile.M);
                nMi = length(M); Midx = NaN(nMi);
                for j = 1:nMi
                    for n = 1:nM
                        if M(j) == obj.meanprofile.M(n)
                            Midx(j) = n;
                            if M(j)<0 || M(j) > 0.9
                                error('Invalid or out-of-bounds Mach number');
                            end
                            break
                        end
                    end
                end
                if isnan(Midx(nMi))
                    error(" Nessun dato trovato per il mach imposto")
                end
            else
                M = obj.meanprofile.M;
                Midx = 1:length(M);
            end
            if nargin == 5
                Fcalc = 0; % se non viene passato Fcalc si suppone che non si voglia forzare il calcolo dei coefficienti
            end
            if  isequal('Landing',flg) || isequal('Take-Off',flg)
                flg = 'High Lift';
            end

            switch flg
                case 'clean'
                    profClass = tridProfInit(obj); % Inizializza l'oggetto profClass con i dati del profilo medio
                    profClass = obj.from2Dto3D(profClass,M,Midx); % Assegna i valori dell'ala 3D
                    profClass.deltaS = 0; profClass.deltaF = 0;

                case 'High Lift'
                    % Include caso di Take-Off e Landing

                    % Cerca se esistono già le informazioni sull'ala pulita
                    % a dato mach
                    calc = 0;
                    if isempty(obj.wing3Ddata) || Fcalc == 1
                        % wing3Ddata è vuoto quindi sicuramente non
                        % conterrà i dati dell'ala pulita
                        calc = 0;
                        nL = 1; % se wing3Ddata è vuoto deve creare come primo elemento l'ala pulita
                    else
                        nL = length(obj.wing3Ddata);
                        for i=1:nL
                            % Controlla se sono già state calcolate le grandezze per l'ala 3D pulita:
                            % preliminarmente è stato già controllat che esistano le grandezze analoghe
                            % 2D e, per come è organizzata la classe, queste occuperanno sempre la
                            % stessa posizione data da Midx
                            if obj.wing3Ddata(i).M(Midx) == M
                                calc = i;
                                oidx = i;
                                %nL = nL + 1;
                                break
                            end
                        end
                    end
                    % Se non sono state trovati i dati dell'ala pulita a
                    % dato Mach li calcola
                    if calc == 0
                        if nL > 1
                            % Se wing3Ddata non è vuoto ma non contiene le
                            % informazioni al mach richiesto allora si crea
                            % un nuovo elemento dell'array dopo l'ultimo
                            % elemento
                            nL = nL + 1;
                        end
                        % Usiamo nL perchè potrebbero esserci già altri
                        % risultati a mach diversi
                        obj.wing3Ddata(nL) = tridProfInit(obj); %Inizializza
                        obj.wing3Ddata(nL) = obj.from2Dto3D(obj.wing3Ddata(nL),M,Midx);
                        % Le informazioni del Profilo con ipersostentatori
                        % sarà salvato in obj.wing3Ddata(nL+1)
                        oidx = nL;
                    end
                    profClass = tridProfInit(obj); % inizializza
                    profClass = obj.from2Dto3D(profClass,M,Midx);
                    % Calcolo degli Effetti dei Flaps

                    if deltaF > 0
                        [cfavg,dcl0_mean,dcl0_tot,cbaroc_avg,cfoc_avg,dclmax_tot] = ...
                            obj.flapEffects(deltaF,Midx,oidx); %Calcola i coefficienti 3D con i flaps
                    else
                        [cfavg,~,~,cbaroc_avg,cfoc_avg,~] = ...
                            obj.flapEffects(deltaF,Midx,oidx); %Calcola i coefficienti 3D con i flaps
                        dcl0_mean = 0; dclmax_tot = 0;
                    end
                    obj.cfocAvg           = cfavg;
                    profClass.a(Midx)     = obj.wing3Ddata(oidx).a(Midx)*( 1 + dcl0_tot/dcl0_mean*...
                        ( cbaroc_avg* (1-cfoc_avg*sin(deltaF*pi/180).^2) -1 ) );
                    profClass.cl0(Midx)   = obj.wing3Ddata(oidx).cl0(Midx) + dcl0_tot;
                    profClass.clmax(Midx) = obj.wing3Ddata(oidx).clmax(Midx) + dclmax_tot;

                    %Calcolo Effetto degli Slat
                    if deltaS > 0
                        [csavg,dCLmaxTotSlat,dClmax2DSlats] = slatEffects(obj,deltaS);
                    else
                        [csavg,~,~] = slatEffects(obj,deltaS);
                        dCLmaxTotSlat = 0;
                        dClmax2DSlats = 0;
                    end
                    obj.csocAvg = csavg;
                    profClass.clmax(Midx) = profClass.clmax(Midx) + dCLmaxTotSlat;

                    profClass.alphamax(Midx)  = ...
                        ( profClass.clmax(Midx) - profClass.cl0(Midx) )/profClass.a(Midx) + ...
                         dAlphaMaxFun(obj.panels(1).sweep,obj.meanprofile.dY) * (deltaS>0);
                    profClass.alphastar(Midx) = profClass.alphamax(Midx)...
                        - ( obj.wing3Ddata(oidx).alphamax(Midx) - obj.wing3Ddata(oidx).alphastar(Midx) );
                    profClass.clstar(Midx)    = profClass.cl0(Midx) + profClass.a(Midx) * profClass.alphastar(Midx);
                    profClass.alpha0l(Midx)   = -profClass.cl0(Midx)/profClass.a(Midx);
                    %
                    SflapTot = 0; SslatTot = 0; cbarocSlat = 0;
                    for i=1:obj.nflaps
                        SflapTot = SflapTot + obj.flaps(i).S;
                        SslatTot = SslatTot + obj.slats(i).S;
                        cbarocSlat = cbarocSlat + ...
                            (obj.slats(i).root.cextoc + obj.slats(i).tip.cextoc)/(obj.nslats*2);
                    end
                    obj.Sflaps = SflapTot; obj.Sslats = SslatTot;
                    dClmax2DSlats = dClmax2DSlats/SslatTot;

                    % Effetti sul Momento
                        % Flaps
                    profClass.cmac(Midx)  = obj.wing3Ddata(oidx).cmac(Midx) + ...
                        Delta_Cm_flaplan_fun(dclmax_tot,cbaroc_avg,cfavg,profClass.clmax(Midx),...
                        obj.Sw, SflapTot,obj.AR,dcl0_mean,...
                        obj.sweepChange(obj.panels(1).sweep,0,0.25,obj.AR,obj.TR),...
                        deltaF,2*(obj.flaps(end).tip.yglob - obj.flaps(1).root.yglob)/obj.bw,obj.TR);
                        % Slats
                    profClass.cmac(Midx) = profClass.cmac(Midx) + ...
                        Delta_Cm_flaplan_fun(dCLmaxTotSlat,cbarocSlat,csavg,profClass.clmax(Midx),...
                        obj.Sw, SslatTot,obj.AR,dClmax2DSlats,...
                        obj.sweepChange(obj.panels(1).sweep,0,0.25,obj.AR,obj.TR),...
                        deltaS,2*(obj.slats(end).tip.yglob - obj.slats(1).root.yglob)/obj.bw,obj.TR);
                    % Rendiamo NaN i valori ottenuti per Mach diversi da
                    nEl = length(profClass.a);
                    idx = 1:nEl; idx(idx == Midx ) = 0;
                    for i = 1:nEl
                        profClass.M(idx>0) = NaN; profClass.a(idx>0) = NaN;
                        profClass.cl0(idx>0) = NaN; profClass.clstar(idx>0) = NaN; 
                        profClass.clmax(idx>0) = NaN; profClass.alphamax(idx>0)= NaN; 
                        profClass.alphastar(idx>0)= NaN;  profClass.alpha0l(idx>0)= NaN; 
                        profClass.cmac(idx>0)= NaN; 
                    end
                    % Assegnamo nome flag
                    if deltaF>0 && deltaS >0
                        profClass.flag = 'Wing with Flaps and Slats';
                        profClass.deltaS = deltaS; profClass.deltaF = deltaF;
                    elseif deltaF>0
                        profClass.flag = 'Wing with Flaps';
                        profClass.deltaS = 0; profClass.deltaF = deltaF;
                    else
                        profClass.flag = 'Wing with Slats';
                        profClass.deltaS = deltaS; profClass.deltaF = 0;
                    end
            end
        end

        function obProf = from2Dto3D(obj,obProf,M,Midx)
            %from2Dto3D: funzione che dato il numero di Mach di volo
            %calcola i coefficienti aerodinamici 3D dell'ala partendo dalel
            %caratteristiche 2D definite nell'oggetto della classe Wing che
            %si ta chiamando.
            %   obProf: oggetto classe profile nel quale immagazzinare i
            %   dati dell'ala 3D
            %   M: numero di Mach di volo
            %   Midx: vettore di indici che specificano a quale riga
            %   corrispondono i valori per un dato mach (ad esempio M =
            %   [0.1,0.2] Midx = [2,1] vuol dire che i valori per M  =0.1
            %   si troveranno nelal seconda riga di meanprofile.a

            % PER ORA CONSIDERA LO SWEEP DEL PRIMO PANNELLO SOLO

            %Cl_alfa 3D: a, sweep_le, AR, sweep_c2, M
            nM = length(M);
            if nargin == 3
                Midx = 1:nM;
            end

            for n = 1:nM
                obProf.a(Midx(n)) = CL_Alfa_fun(obj.meanprofile.a(Midx(n)), obj.panels(1).sweep,obj.AR,...
                    obj.sweepChange(obj.panels(1).sweep,0,0.5,obj.AR,obj.TR),...
                    M(n));
                %Cl_max 3D: Mean_Cl_max,sweep,dy,c, M_inf: dCLMaxFun
                obProf.clmax(Midx(n)) = CL_max_fun(...
                    obj.meanprofile.clmax(Midx(n)),obj.panels(1).sweep,obj.meanprofile.dY,obj.meanprofile.c,M(n)); % dCLMaxFun ha problemi nell'estrapolare

                obProf.clstar(Midx(n)) = obj.meanprofile.clstar(Midx(n));
                obProf.alpha0l(Midx(n)) = obj.meanprofile.alpha0l(Midx(n));
                obProf.cl0(Midx(n)) = obProf.a(Midx(n)) * (-obProf.alpha0l(Midx(n)));

                %Alfa_max 3D: clmax,a,alpha0l,dY,sweep
                obProf.alphamax(Midx(n)) = AlphaMaxFun(obProf.clmax(Midx(n)),obProf.a(Midx(n)),...
                    obProf.alpha0l(Midx(n)),obj.meanprofile.dY,obj.panels(1).sweep);
                
                %Alfa_star 3D:
                obProf.alphastar(Midx(n)) = ...
                    (obProf.clstar(Midx(n)) - obProf.cl0(Midx(n)))/obProf.a(Midx(n));

            end
            obProf.flag = 'Wing in Clean Configuration';
        end

        function [cfavg,dcl0_mean,dcl0_tot,cbaroc_avg,cfoc_avg,dclmax_tot] = flapEffects(obj,deltaF,Midx,oidx)
            %flapEffects; calcola gli effetti dei flaps sulle
            %caratteristiche 3D dell'ala
            %   deltaF: deflessioen dei flaps in deg
            %   Midx: indice che contiene la posizione delle grandezze
            %       aerodinamiche nell' oggetto profilo che contiene i dati
            %       dell'ala pulita.
            %   oidx: indice che indica a che posizione dell'array
            %       wing3Ddata si trovano i dati dell'ala pulita

            % Inizializzazione delle Variabili di Output
            dcl0_mean = 0; dcl0_tot = 0; cbaroc_avg = 0;
            cfoc_avg = 0; dclmax_tot = 0; weiS = 0; cfavg = 0;
            for i =1:obj.nflaps
                % variabili usate per il calcolo delle
                % caratteristiche high lift
                % (1) dCl0 flap
                [dclo,cboc,alphadf] = dCl02DHLFun(deltaF,obj.flaps(i),'root',Midx);
                obj.flaps(i).root.HLauxVariables('dCl02D', dclo);
                obj.flaps(i).root.HLauxVariables('cbaroc', cboc);
                obj.flaps(i).root.HLauxVariables('alphaDeltaf', alphadf);
                obj.flaps(i).root.HLauxVariables('dClmax2D', dClMaxHLFun(deltaF,...
                    obj.flaps(i),'root'));
                [dclo,cboc,alphadf] = dCl02DHLFun(deltaF,obj.flaps(i),'tip',Midx);
                obj.flaps(i).tip.HLauxVariables('dCl02D', dclo);
                obj.flaps(i).tip.HLauxVariables('cbaroc', cboc);
                obj.flaps(i).tip.HLauxVariables('alphaDeltaf', alphadf);
                obj.flaps(i).tip.HLauxVariables('dClmax2D', dClMaxHLFun(deltaF,...
                    obj.flaps(i),'tip'));

                % Calcolo Pesi Corde Flaps
                k1 = 0.5*obj.flaps(i).root.c*obj.flaps(i).b/obj.flaps(i).S;
                k2 = 0.5*obj.flaps(i).tip.c*obj.flaps(i).b/obj.flaps(i).S;

                % Calcoli effetto flaps su cm
                weiS = k1 + k2 + weiS; % Somma dei pesi
                cfavg = cfavg + k1*obj.flaps(i).root.cfoc + k2*obj.flaps(i).tip.cfoc;
                %dClMax2D_flap
                %dClMax2D medio
                obj.flaps(i).HLauxVariables('dCLmax_mean', k1*obj.flaps(i).root.HLauxVariables('dClmax2D')+...
                    k2*obj.flaps(i).tip.HLauxVariables('dClmax2D'));

                sweep025 = ...
                    obj.sweepChange(obj.panels(1).sweep,0,0.25,obj.AR,obj.TR);
                Ksweep = (1-0.08*cos(sweep025*pi/180).^2)*cos(sweep025*pi/180).^(3/4) ;
                SfoS = 2*obj.flaps(i).S/obj.Sw;
                %dCLmax
                obj.flaps(i).coeffsShiftFun('dCLmax',SfoS*Ksweep*...
                    obj.flaps(i).HLauxVariables('dCLmax_mean'));

                %%dCl02D_flap
                %dCl02D medio
                obj.flaps(i).HLauxVariables('dCl02D_mean', ...
                    k1*obj.flaps(i).root.HLauxVariables('dCl02D')+...
                    k2*obj.flaps(i).tip.HLauxVariables('dCl02D'));
                % a medio flap
                obj.flaps(i).HLauxVariables('a_mean',k1*obj.flaps(i).root.a(Midx)+...
                    k2*obj.flaps(i).tip.a(Midx));
                % alphaDeltaf medio
                obj.flaps(i).HLauxVariables('alphaDeltaf' , ...
                    k1*obj.flaps(i).root.HLauxVariables('alphaDeltaf')+...
                    k2*obj.flaps(i).tip.HLauxVariables('alphaDeltaf'));
                %Kb flap
                obj.flaps(i).HLauxVariables('Kb', kbFun(...
                    [2*obj.flaps(i).root.yglob/obj.bw,...
                    2*obj.flaps(i).tip.yglob/obj.bw],obj.TR));
                %Kc flap
                obj.flaps(i).HLauxVariables('Kc', kcFun( ...
                    obj.flaps(i).HLauxVariables('alphaDeltaf')...
                    ,obj.AR,1));
                %dCL0 flaps
                obj.flaps(i).coeffsShiftFun('dCL0',...
                    obj.wing3Ddata(oidx).a(Midx)/obj.flaps(i).HLauxVariables('a_mean')*...
                    obj.flaps(i).HLauxVariables('dCl02D_mean')*...
                    obj.flaps(i).HLauxVariables('Kb')*...
                    obj.flaps(i).HLauxVariables('Kc') );

                %a 3D con flap
                dcl0_mean = dcl0_mean + ...
                    obj.flaps(i).HLauxVariables('dCl02D_mean')*2*obj.flaps(i).S/obj.Sw;
                dcl0_tot = dcl0_tot + obj.flaps(i).coeffsShiftFun('dCL0');

                cbaroc_avg = cbaroc_avg + (obj.flaps(i).root.HLauxVariables('cbaroc')...
                    + obj.flaps(i).tip.HLauxVariables('cbaroc'))/(2*obj.nflaps);

                cfoc_avg = cfoc_avg + ...
                    (obj.flaps(i).root.cfoc/obj.flaps(i).root.HLauxVariables('cbaroc') +...
                    obj.flaps(i).tip.cfoc/obj.flaps(i).tip.HLauxVariables('cbaroc'))/(2*obj.nflaps);

                dclmax_tot = dclmax_tot + obj.flaps(i).coeffsShiftFun('dCLmax');

            end
            cfavg = cfavg/weiS; % Somma dei pesi
        end

        function [csavg,dCLmaxTotSlat, dClmax2DSlats] = slatEffects(obj,deltaS)
            for i=1:obj.nslats
                % Sezione Slat di Radice
                obj.slats(i).root.HLauxVariables('etaMax',...
                    etaMaxFun(obj.slats(i).root.LERc/obj.slats(i).root.tc,1)); % DatoLER/t potrebbe essere sbagliato
                obj.slats(i).root.HLauxVariables('etaDelta',...
                    etaDeltaFun(deltaS,1));
                obj.slats(i).root.HLauxVariables('CloDs',...
                    clOdsSlat(obj.slats(i).root.cfoc));
                obj.slats(i).root.HLauxVariables('dClmax_slat',...
                    obj.slats(i).root.HLauxVariables('etaMax')*...
                    obj.slats(i).root.HLauxVariables('etaDelta')*...
                    obj.slats(i).root.HLauxVariables('CloDs')*...
                    deltaS*obj.slats(i).root.cextoc);
                % Sezione Slat di Estremità
                obj.slats(i).tip.HLauxVariables('etaMax',...
                    etaMaxFun(obj.slats(i).tip.LERc/obj.slats(i).tip.tc,1));
                obj.slats(i).tip.HLauxVariables('etaDelta',...
                    etaDeltaFun(deltaS,1));
                obj.slats(i).tip.HLauxVariables('CloDs',...
                    clOdsSlat(obj.slats(i).tip.cfoc));
                obj.slats(i).tip.HLauxVariables('dClmax_slat',...
                    obj.slats(i).tip.HLauxVariables('etaMax')*...
                    obj.slats(i).tip.HLauxVariables('etaDelta')*...
                    obj.slats(i).tip.HLauxVariables('CloDs')*...
                    deltaS*obj.slats(i).tip.cextoc);
            end

            % Calcolo degli effetti 3D
            dCLmaxTotSlat = 0; weiS = 0; csavg = 0; dClmax2DSlats = 0;
            for i=1:obj.nslats
                k1 = 0.5*obj.slats(i).root.c*obj.slats(i).b/obj.slats(i).S;
                k2 = 0.5*obj.slats(i).tip.c*obj.slats(i).b/obj.slats(i).S;

                weiS = k1 + k2 + weiS; % Somma dei pesi
                csavg = csavg + k1*obj.slats(i).root.cfoc + k2*obj.slats(i).tip.cfoc;

                obj.slats(i).HLauxVariables('dCmax2D_slat',...
                    k1*obj.slats(i).root.HLauxVariables('dClmax_slat') + ...
                    k2*obj.slats(i).tip.HLauxVariables('dClmax_slat') );
                sweep025 = ...
                    obj.sweepChange(obj.panels(1).sweep,0,0.25,obj.AR,obj.TR);
                Ksweep = (1-0.08*cos(sweep025*pi/180).^2)*cos(sweep025*pi/180).^(3/4) ;
                obj.slats(i).coeffsShiftFun('dCLmax_slat',...
                    obj.slats(i).HLauxVariables('dCmax2D_slat')*Ksweep*2*obj.slats(i).S/obj.Sw );
                dCLmaxTotSlat = dCLmaxTotSlat + obj.slats(i).coeffsShiftFun('dCLmax_slat');

                dClmax2DSlats = dClmax2DSlats + obj.slats(i).HLauxVariables('dCmax2D_slat')*obj.slats(i).S;

            end
            csavg = csavg/weiS;
        end

        function [eps0,depsoda] = downwashEstimation(~,Xaw,Zaw,Xah,Zah,iw,aol,b,AR,a,sweep)
            %downwashEstimation: valuta il gradiente di downwash agente sul
            %piano di coda
            % Xaw, Zaw:     posizione del centro aerodinamico dell'ala
            % Xah, Zah:     posizione del centro aerodinamico del piano
            %                   orizzontale
            % iw:           angolo calettamento dell'ala
            % aol:          angolo portanza nulla dell'ala
            % b:            apertura alare
            % AR:           aspect ratio dell'ala
            % a:            pendenza retta di portanza dell'ala in [deg^-1]
            % sweep:        Angolo di sweep dell'ala in deg

            sweep = sweep*pi/180; % Porta lo sweep in radianti
            d = sqrt((Xah - Xaw)^2 + (Zah - Zaw)^2);
            phi = atan((Zah - Zaw)/(Xah - Xaw))*180/pi + iw - aol;
            r = Xah - Xaw; r = 2*r/b;
            m = d*sin(phi*pi/180); m = 2*m/b;
            % Coefficiente di Downwash senza freccia
            Keps = r/(r^2+m^2) * 0.4876/sqrt(r^2 + 0.6319+m^2) + ...
                ( 1 + (r^2 / (r^2 + 0.79 + 5.0734*m^2) )^(0.3113) )*...
                (1 - sqrt( m^2/(1 + m^2) ) );
            % Correzione a causa della freccia dell'ala
            KswoKeps = ( (0.1124 + 0.1265*sweep + 0.1766*sweep^2)/r^2 + 0.1124/r + 2 )/...
                (0.1124/r^2 + 0.1124/r + 2);
            % Coefficiente di Downwash corretto
            Keps = Keps*KswoKeps;

            depsoda = Keps*(a*180/pi)/(pi*AR);
            eps0 = depsoda*aol;
        end

        %% Funzioni di Drag Estimation

        function [Cd0,dCd0,e] = DragEstLifSur(obj,profObj,fusobj,h,Neng,QS,QSN,M)
            %DragEstLifSur calcola il Cd0 e il fattore di Oswald per una
            %superficie portante
            %   fusobj: oggetto di tipo fusoliera i cui elementi sono i
            %       DIAMETRI di fusoliera
            %   profObj: oggetto della classe profilo, nel quale sono
            %       contenuti i dati aerodinamici a mach fissato
            %   h: quota di volo
            %   fCond: flag che indica la cond. di volo ( 'Cruise'
            %       'TakeOff' 'Landing')
            %   Neng: Numero di Motori
            obj.tau = obj.panels(end).tip.tc/obj.panels(1).root.tc;
            profObj.K = obj.kdragWing;
            profObj.h = h(:);
            sweepXT = obj.sweepChange(obj.panels(1).sweep,0,...
                obj.meanprofile.xtc,obj.AR,obj.TR);
            % Calcolo della corda esposta
            %fusobj.dfmax
            for i =1:obj.npanels
                crootexp = obj.panels(i).panelInterp(fusobj.dfmax*0.5,1); % sarebbe il raggio massimo di fusoliera
                if ~isempty(crootexp)
                    npan = i;
                    break
                end
            end
            % Calcolo dell'area bagnata
            Sexl = 0;
            for j = 1:npan
                Sexl = Sexl + obj.panels(j).S; % Area dei pannelli coperti dalla fusoliera
            end
            % Calcolo della Superficie Esposta: Si sottrae alla superficie
            % dell'ala quella dei pannelli che sono anche parzialmente
            % inglobat nella fusoliera e si aggiunge la frazione di area
            % del pannello che è parzialmente inglobato
            Sexp = obj.Sw - 2*( Sexl - 0.5*(obj.panels(npan).tip.yglob - fusobj.dfmax*0.5)*(crootexp+obj.panels(npan).tip.c) );
            % Calcolo Area Bagnata
            obj.Swet = 2*Sexp*( 1 + 0.25*obj.panels(1).root.tc*( 1 + obj.tau*obj.TR )/( 1 + obj.TR) );

            if nargin == 7
                % Se non viene dato il Mach calcola con tutti i profili a
                % disposizione
                M = profObj.M;
                nM = length(M);
                Midx = 1:nM;
            else
                nM = length(profObj.M);
                for j = 1:nM
                    if profObj.M(j) == M
                        Midx = j;
                        nM = 1;
                        break;
                    end
                end
            end
            e     = zeros(nM,1);
            Cd0   = zeros(nM,1);
            dCd0  = zeros(nM,1);
            for n=1:nM
                [~, a, ~, rho,~,mu] = atmosisa(h(n)); % ATTENZIONE: restituisce mu solo con MATLAB 2024b
                profObj.Re(Midx(n)) = rho*a*profObj.M(Midx(n))*obj.meanprofile.c/mu;
                if profObj.M(Midx(n))<0.9
                    ReCutOff = 38.21*( obj.meanprofile.c/profObj.K)^1.053;
                else
                    ReCutOff = 44.62*(( obj.meanprofile.c/profObj.K)^1.053 )*profObj.M(Midx(n))^1.16;
                end
                CfLam = 1.328/sqrt(profObj.Re(Midx(n))); %Cf laminare da Blasius
                if profObj.Re(Midx(n)) < ReCutOff % Reynolds per il quale il Cf è costante nell'abaco di Moody?
                    CfTurb = 0.455/( (log10(profObj.Re(Midx(n)))^2.58) * (1+0.144*profObj.M(Midx(n))^2)^0.65 );
                else
                    CfTurb = 0.455/( (log10(ReCutOff)^2.58) * (1+0.144*profObj.M(Midx(n))^2)^0.65 );
                end
                CfUp = obj.meanprofile.xtrUp*CfLam + (1 - obj.meanprofile.xtrUp)*CfTurb;
                CfLow = obj.meanprofile.xtrLow*CfLam + (1 - obj.meanprofile.xtrLow)*CfTurb;
                profObj.Cf(Midx(n)) = CfUp*0.5 + CfLow*0.5;
                profObj.FF(Midx(n)) = ( 1 + 0.6/obj.meanprofile.xtc * obj.meanprofile.tc + 100*obj.meanprofile.tc^4 ) * ...
                    ( 1.34*profObj.M(Midx(n))^0.18 * cos(sweepXT*pi/180)^0.28 ); % Form Factor
                profObj.Cd0(Midx(n)) = obj.Swet/obj.Sw * profObj.Cf(Midx(n)) * profObj.FF(Midx(n)) * QS * QSN; % 1, 1.3 cono fattori di interf. da definire

                if profObj.deltaF ~= 0
                    dCd0(Midx(n)) = obj.HLDrag(obj.cfocAvg, obj.Sflaps*2,profObj.deltaF);
                else
                    dCd0(Midx(n)) = 0;
                end
                if profObj.deltaS ~= 0
                    dCd0(Midx(n)) = dCd0(Midx(n)) + obj.HLDrag(obj.csocAvg, obj.Sslats*2,profObj.deltaS);
                end
                fl = 0.005*( 1+1.5*(obj.TR-0.6)^2 );
                AReff = (obj.bw - fusobj.dfmax)^2 / Sexp;
                profObj.e(Midx(n)) = 1/( (1+0.12*profObj.M(Midx(n))^6) *(1 + (0.142+fl*AReff*(10*obj.meanprofile.tc)^0.33)/...
                    (cos(obj.sweepChange(obj.panels(1).sweep,0,0.25,obj.Sw,obj.TR)*pi/180))^2 + ...
                    0.1*(Neng*3+1)/(4+AReff)^0.8) );
                %if obj.bw ~= 0
                % Effetto delle Winglet
                profObj.e(Midx(n)) = profObj.e(Midx(n))*(1+2*obj.hwglt/obj.bw)^2;
                %end
                e(Midx(n)) = profObj.e(Midx(n));
                Cd0(Midx(n)) = profObj.Cd0(Midx(n));
                % profObj.dCd0(Midx(n)) = dCd0;
            end
        end

        function dC0 = HLDrag(obj,cfoc,ShL,d)
            %HLDrag: calcola l'aumento di Cd dovouto agli ipersostentatori
            %    cfoc: valore medio della corda dell' HL sulla corda del profilo
            %    ShL: Valore totale dell'area dei flap

            dC0 = 0.0074*cfoc*(ShL/obj.Sw)*(d - 10);
        end

        function K = kdragWing(obj)
            K = 0.00635*1e-3;% m
            %             switch obj.liftSurfFlag
            %                 case 'Wing'
            %
            %                 case 'Horizontal'
            %
            %                 case 'Vertical'
            %
            %             end
        end

    end
end
