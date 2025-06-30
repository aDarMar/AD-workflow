classdef WingClass %< handle %<WingClass è una sottoclassed della classe predefinita handle
    properties
        %% Geometry
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
        %% Aerodynamic Data
        wing3Ddata
        prf3DClean
        % Dati Aggiuntivi Piano di Coda
        eps0 {mustBeNumeric}    %angolo di downwash per alpha = 0
        depsda {mustBeNumeric}  %Gradiente di downwash
        % High Lift Devices
        prf3Dflaps
        prf3Dflapslat

        %% Global Coordinates 
        Xapex {mustBeNumeric}
        Yapex {mustBeNumeric}
        Zapex {mustBeNumeric}
        TST
    end
    methods
        function obj = WingClass(bs,sweeps,dihedrals,iang,apexC,M,...
                sectsGeom, sectsAero,...
                HLflag,sectsHL,...
                misc,tFL ) %Costruttore
            %INPUT
            %   HLflag: variable containing the type of HL surface, for
            %       example if there are three flaps, it will contain three
            %       chars 'flap'
            %   sectsHL: for flaps -> [ 2*y/b,cf/c,flap_type ]
            %            for slats -> [ 2*y/b, cs/c, c_ext/c ]
            %   tFL: a flag that if it passed tells the program to consider
            %       the extended aerodynamics inputs
            if nargin == 1
                % Lettura XML
            else
                obj.TST = 1;
                nM = length(M); % Numero di condizioni di volo
%                 % Assegnazione Coordinate Apice Superficie
%                 obj.Xapex = apexC(1);
%                 obj.Yapex = apexC(2);
%                 obj.Zapex = apexC(3);
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
            % Assign rhe coordinates of the wing apex to the innermost
            % profile
            i = 1;
            obj.panels(i).root.xglob = apexC(1);
            obj.panels(i).root.yglob = apexC(2);
            obj.panels(i).root.zglob = apexC(3);
            
            [obj.panels(i).tip.xglob,obj.panels(i).tip.yglob,obj.panels(i).tip.zglob] ...
                = obj.panels(i).globalCoords( obj.panels(i).root.xglob,...
                    obj.panels(i).root.yglob,obj.panels(i).root.zglob );  
                
            for i=2:obj.npanels
                % Assign to outer's root the values of inner's tip
                obj.panels(i).root.xglob = obj.panels(i-1).tip.xglob;
                obj.panels(i).root.yglob = obj.panels(i-1).tip.yglob;
                obj.panels(i).root.zglob = obj.panels(i-1).tip.zglob;
                % Calculate the values of tip
                [obj.panels(i).tip.xglob,obj.panels(i).tip.yglob,obj.panels(i).tip.zglob] ...
                    = obj.panels(i).globalCoords( obj.panels(i-1).tip.xglob,...
                    obj.panels(i-1).tip.yglob,obj.panels(i-1).tip.zglob );
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
                obj.wing3Ddata  = ProfileClass.empty;
                obj.flaps       = PanelClass.empty;
                obj.slats       = PanelClass.empty;
                obj.commandSurf = PanelClass.empty;
                sectsHL         = [sectsHL(:,1)*0.5*obj.bw,sectsHL]; % Trasforma le coordinate adimensionali in dimensionali
                m = 1; mm = 1;
                while mm < length(HLflag) + 1
                    switch HLflag{mm}
                        case 'flaps'
                            obj.nflaps = obj.nflaps + 1;
                            % Dati Geometrici delle sezioni ei profili
                            varG = NaN(2,8+3);
                            % 2 coordinate di flap e 8 valori
                            varG(:,9:11) = sectsHL(m:m+1,2:end);   % copia le grandezze di input dei flaps nella var temp
                            varG(:,9)    = varG(:,9)*0.5*obj.bw;   % trasforma dy/b in coord. dimensionali
                            % Dati Aerodinamici delle Sezioni dei Profili
                            if nargin < 12 || tFL == 0 % this means that tFL has not been defined  
                                [varG,varA,coords]    = obj.HLAssign( varG,nM );
                            else
                                [varG,varA,coords]    = obj.HLAssign( varG,nM,9 );
                            end
                            obj.flaps(obj.nflaps) = PanelClass(varG(2,9) - varG(1,9),0,0,M,...
                                varG(1,1:8),varA(1:nM,:),varG(2,1:8),varA(nM+1:2*nM,:),HLflag{mm},sectsHL(m:m+1,:));
                            % Assegna le coordinate globali
                            obj.flaps(obj.nflaps).root.xglob = coords(1,1); obj.flaps(obj.nflaps).root.zglob = coords(1,2);
                            obj.flaps(obj.nflaps).tip.xglob  = coords(2,1); obj.flaps(obj.nflaps).tip.zglob = coords(2,2);
                            m = m + 2; mm = mm + 1;
                        case 'slats'
                            obj.nslats = obj.nslats + 1;
                            % Dati Geometrici delle sezioni ei profili
                            varG = NaN(2,8+3);
                            % 2 coordinate di flap e 8 valori
                            varG(:,9:11) = sectsHL(m:m+1,2:end);   % copia le grandezze di input dei flaps nella var temp
                            varG(:,9)    = varG(:,9)*0.5*obj.bw;   % trasforma dy/b in coord. dimensionali
                            % Dati Aerodinamici delle Sezioni dei Profili
                            if nargin < 12 || tFL == 0 % this means that tFL has not been defined  
                                [varG,varA,coords]    = obj.HLAssign( varG,nM );
                            else
                                [varG,varA,coords]    = obj.HLAssign( varG,nM,9 );
                            end
                            obj.slats(obj.nslats) = PanelClass(varG(2,9) - varG(1,9),0,0,M,...
                                varG(1,1:8),varA(1:nM,:),varG(2,1:8),varA(nM+1:2*nM,:),HLflag{mm},sectsHL(m:m+1,:));
                            % Assegna le coordinate globali
                            obj.slats(obj.nslats).root.xglob = coords(1,1); obj.slats(obj.nslats).root.zglob = coords(1,2);
                            obj.slats(obj.nslats).tip.xglob  = coords(2,1); obj.slats(obj.nslats).tip.zglob = coords(2,2);
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
            % Specifies which method shall be used to evaluate Alpha0L and
            % CM of the mean profile
            if length( sectsAero(1,:) ) > 8
                met_flg = 2;
            else
                met_flg = 1;
            end
            [tmp1,tmp2]     = obj.meanProfileMod(met_flg);
            obj.meanprofile           = ProfileClass(tmp1,tmp2,M);
            [~,xm,ym,zm]              = obj.macCalc;
            obj.meanprofile.xglob     = xm;
            obj.meanprofile.yglob     = ym;
            obj.meanprofile.zglob     = zm;


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
                mac  = mac + obj.panels(i).mac * 2*obj.panels(i).S/obj.Sw; %il 2 sta perchè panel.S è l'area del pannello e wing.S è l'area totale dell'ala
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
        
        function [vout,vout2] = meanProfileMod(obj,method)
            %meanProfile: calcola le grandezze medie dell'ala.
            % Funzione che calcola il profilo medio. It can employ two
            % methods:
            %   1. Weightned averanges with chords
            %   2. Trapezoidal Integration
            %INPUT
            %   method: flag that specifies the method for obtaining the
            %       alpha0L and CM
            %OUTPUT:
            %   vout2: aerodynamic data of mean profile. Angles and
            %       derivatives are in [deg] and [1/deg] respectively
            av  = NaN(1,obj.npanels+1);
            dYv = av;xrtcv = av; xtrUpv = av; xrtLowv = av; tcv = av;
            nM  = length(obj.panels(1).root.a);
            av  = NaN(nM,obj.npanels+1);
            alphamaxv = av; alpha0lv = av; cmacv = av; alphastarv = av;
            cl0v = av; clstarv = av; clmaxv = av;
            
            for i = 1:obj.npanels
                tcv(i)         = obj.panels(i).root.tc;
                xrtcv(i)       = obj.panels(i).root.xtc;
                xtrUpv(i)      = obj.panels(i).root.xtrUp;
                xrtLowv(i)     = obj.panels(i).root.xtrLow;
                dYv(i)         = obj.panels(i).root.dY;
                
                av(:,i)        = obj.panels(i).root.a(:);
                cl0v(:,i)      = obj.panels(i).root.cl0(:);
                clstarv(:,i)   = obj.panels(i).root.clstar(:);
                clmaxv(:,i)    = obj.panels(i).root.clmax(:);
                alphamaxv(:,i) = obj.panels(i).root.alphamax - obj.panels(i).root.eps(:);
                %alpha0lv(:,i) = obj.panels(i).root.alpha0l - obj.panels(i).root.eps(:);
                alphastarv(:,i)= (obj.panels(i).root.clstar(:) -  obj.panels(i).root.cl0(:))./obj.panels(i).root.a(:);
                %cmacv(:,i) = obj.panels(i).root.cmac(:);
            end
            i = obj.npanels;
            tcv(i+1)     = obj.panels(i).tip.tc;
            xrtcv(i+1)   = obj.panels(i).tip.xtc;
            xtrUpv(i+1)  = obj.panels(i).tip.xtrUp;
            xrtLowv(i+1) = obj.panels(i).tip.xtrLow;
            dYv(i+1)     = obj.panels(i).tip.dY;
            
            av(:,i+1) = obj.panels(i).tip.a;
            cl0v(:,i+1) = obj.panels(i).tip.cl0;
            clstarv(:,i+1) = obj.panels(i).tip.clstar;
            clmaxv(:,i+1) = obj.panels(i).tip.clmax;
            alphamaxv(:,i+1) = obj.panels(i).tip.alphamax - obj.panels(i).tip.eps;
            %alpha0lv(:,i+1) = obj.panels(i).tip.alpha0l - obj.panels(i).tip.eps;
            alphastarv(:,i+1) = (obj.panels(i).tip.clstar -  obj.panels(i).tip.cl0)./obj.panels(i).tip.a;
            %cmacv(:,i+1) = obj.panels(i).tip.cmac;
            [cm,xLE_MAC,~,~] = macCalc(obj);
            
            
            if nargin <2 || method == 1
                % Weightned averanges.
                for i = 1:obj.npanels
                    alpha0lv(:,i) = obj.panels(i).root.alpha0l - obj.panels(i).root.eps(:);
                    cmacv(:,i)    = obj.panels(i).root.cmac(:);
                end
                i = obj.npanels;
                alpha0lv(:,i+1) = obj.panels(i).tip.alpha0l - obj.panels(i).tip.eps;
                cmacv(:,i+1)    = obj.panels(i).tip.cmac;
                vout2 = [obj.weightAvg(av),...
                    obj.weightAvg(cl0v),obj.weightAvg(clstarv),obj.weightAvg(clmaxv),...
                    obj.weightAvg(alphamaxv),obj.weightAvg(alpha0lv),...
                    obj.weightAvg(alphastarv),obj.weightAvg(cmacv)];
            else
                % Trapezoidal Integration: changes only alpha_0L and CM
                %% Definition of Spanwise Sections
                dy = 0.5; %dy = floor(0.5*obj.bw/dy);
                yvec  = 0:dy:obj.panels(1).tip.yglob;
                yvec2 = obj.panels(2).root.yglob:dy:obj.panels(2).tip.yglob;
                yvec  = [yvec(1:end-1),yvec2(1:end-1),obj.panels(2).tip.yglob];
                %yvec = [0,1,2,3,4,4.333458599,5,6,7,8,9,10,11,12,13,14,15,16,17,17.33383439];
                n_stats = length(yvec);
                %% Building Interpolation Vectors
                cvet     = interp1( [obj.panels(1).root.yglob,obj.panels(2).root.yglob,obj.panels(2).tip.yglob],...
                    [obj.panels(1).root.c,obj.panels(2).root.c,obj.panels(2).tip.c],yvec );
                eps_vet  = interp1( [obj.panels(1).root.yglob,obj.panels(2).root.yglob,obj.panels(2).tip.yglob],...
                    [obj.panels(1).root.eps,obj.panels(2).root.eps,obj.panels(2).tip.eps],yvec );
                azl_vet  = interp1( [obj.panels(1).root.yglob,obj.panels(2).root.yglob,obj.panels(2).tip.yglob],...
                    [obj.panels(1).root.alpha0l,obj.panels(2).root.alpha0l,obj.panels(2).tip.alpha0l], yvec );
                CMac_vet = interp1( [obj.panels(1).root.yglob,obj.panels(2).root.yglob,obj.panels(2).tip.yglob],...
                    [obj.panels(1).root.cmac,obj.panels(2).root.cmac,obj.panels(2).tip.cmac],yvec );
                Cla_vet  = interp1( [obj.panels(1).root.yglob,obj.panels(2).root.yglob,obj.panels(2).tip.yglob],...
                    [obj.panels(1).root.a,obj.panels(2).root.a,obj.panels(2).tip.a],yvec );
                Xle_vec  = interp1( [obj.panels(1).root.yglob,obj.panels(2).root.yglob,obj.panels(2).tip.yglob],...
                    [obj.panels(1).root.xglob,obj.panels(2).root.xglob,obj.panels(2).tip.xglob],yvec );
                Xac_vet  = interp1( [obj.panels(1).root.yglob,obj.panels(2).root.yglob,obj.panels(2).tip.yglob],...
                    [obj.panels(1).root.x_ac,obj.panels(2).root.x_ac,obj.panels(2).tip.x_ac],yvec );
                %% Wing Alpha-Zero Lift and Aerodynamic Twist     
                azl_int  =  cvet.*( azl_vet-eps_vet );
                eps_int  =  ( -azl_vet(:) + azl_vet(1) ).*Cla_vet(:).*cvet(:);         
                azl_mean = obj.trapezoidal_int( azl_int,yvec ); % Zero-Lift Angle [deg]
                azl_mean = azl_mean*2/obj.Sw;
                eps_a    = obj.trapezoidal_int( eps_int,yvec ); % Mean Aerodynamic Twist. Check also after vout2
                
                %% Mean Aerodynamic Twist
                
                %% Moments and Wing Aerodynamic Center
                % Vector of coordinates of wing aerodynamic centers
                % xc_4 = X_le + xac,p(y)/c(y)* c(y)
                xc_4    =  Xac_vet.*cvet + Xle_vec;       
                % % Distance between wing and profiles aerodynamic center
                alpha   = azl_mean:0.5:obj.weightAvg(alphastarv); % Sweep in alpha from a0l to alpha*
                alpha = [-3,alpha];
                % The x_ac/wing is calculated using the mean airfoil
                % Cl_alpha
                [xc_axw,Cm0,cm_add] = cm_alpha( obj,alpha,obj.weightAvg(av),...
                    yvec,cvet,Xle_vec,Cla_vet,eps_vet,azl_vet,CMac_vet,Xac_vet,cm,xLE_MAC );
                % REMARK: cm_add will be equal to Cma because they are both
                % calculated @wing mac
                % Additional Load@CL = 0
                x_ref   = ( xLE_MAC + cm*xc_axw );
                x1      =  x_ref - xc_4; 
                Cma_int = ( azl_mean+eps_vet-azl_vet ).*Cla_vet.*cvet.*x1;
                Cma     = obj.trapezoidal_int( Cma_int,yvec );
                Cma     = 2*Cma/(obj.Sw*cm);
                % We found Cm@CL = 0 that is the CM@ac_wing: THIS IS TRUE
                % ONLY WHEN THE REFERENCE AXIS IS THE WING MAC

                vout2 = [obj.weightAvg(av),...
                    obj.weightAvg(cl0v),obj.weightAvg(clstarv),obj.weightAvg(clmaxv),...
                    obj.weightAvg(alphamaxv),azl_mean,...
                    obj.weightAvg(alphastarv),Cm0+Cma,xc_axw,0];
                eps_a = 2*eps_a/( vout2(1)*obj.panels(end).tip.c*obj.panels(end).tip.yglob ); % eps_a = sum/( Cla_avg*c_tip*b/2 )
                vout2 = [vout2,eps_a];
            end
            vout = [cm,nan,obj.weightAvg(tcv),nan,obj.weightAvg(xrtcv),obj.weightAvg(xtrUpv),...
                obj.weightAvg(xrtLowv),obj.weightAvg(dYv)];

        end

        function polydrag = poly_drag( obj,alpha,cds )
            %poly_drag: function that calculates the mean cd for each alpha
            % and interpolates Cd - alpha values giving back a polyfit object
            %INPUT:
            %   alpha: column vector of alphas at which the cd are
            %       calculated;
            %   cds: vector containing for each row the values of cd for
            %       root kinik and tip at a given alpha.

            % Regression to find the experimental Cd-alpha values
            n_alpha = length( alpha( : ) );
            cd_av = nan(n_alpha,1);
            if log10( cds(1,1) ) > 0
                % Check to see if Cds are given as drag numbers or
                % naturally
                cds = cds*1e-4; % fromt drag count to normal scale
            end
            for i = 1:n_alpha
                cd_av(i) = obj.weightAvg( cds(i,:) );
            end
            % Curve fitting of drag
            n = 5;
            polydrag = polyfit( alpha,cd_av,n );
        end
        
        function [x_ref,cm0,cm_add] = cm_alpha( obj,alpha,Cla,yvec,cvet,x_vet,cla_vet,eps_vet,azl_vet,cm0_vet,xac_vet,MAC,XleMAC )
            %cm_alpha: function that evaluates the Cm-alpha curve for a
            %given set of alpha using cm data from sections
            %INPUT
            %   alpha: vector of alpha to calculate Cm. They should be in
            %       the lienar range of the sections Cl-alpha
            %   Cla: wing lift slope
            %   yvec ... XleMAC: data associated with sections. If they are
            %       omitted, the program will calculate them by itself
            %OUTPUT
            %   xoc_ref: wing x_ac position as a fraction of MAC
            %   cm_0: pitching moment due to basic load
            %   cm_add: vector of cm due to additional load evaluated at
            %           at alpha given in alpha vector. BE CAREFUL: the
            %           function exits when it finds the x_ac,wing so it is
            %           expected that cm_add vs alpha will be all the same
            if nargin < 3
                % Wing stencil Definition
                dy = 0.5; %dy = floor(0.5*obj.bw/dy);
                yvec  = 0:dy:obj.panels(1).tip.yglob;
                yvec2 = obj.panels(2).root.yglob:dy:obj.panels(2).tip.yglob;
                yvec  = [yvec(1:end-1),yvec2(1:end-1),obj.panels(2).tip.yglob];
                % Interpolation of Needed data
                % chords
                cvet     = interp1( [obj.panels(1).root.yglob,obj.panels(2).root.yglob,obj.panels(2).tip.yglob],...
                    [obj.panels(1).root.c,obj.panels(2).root.c,obj.panels(2).tip.c],yvec );
                % leading edges x coords
                x_vet   = interp1( [obj.panels(1).root.yglob,obj.panels(2).root.yglob,obj.panels(2).tip.yglob],...
                    [obj.panels(1).root.xglob,obj.panels(2).root.xglob,obj.panels(2).tip.xglob],yvec );
                % cl_a
                cla_vet = interp1( [obj.panels(1).root.yglob,obj.panels(2).root.yglob,obj.panels(2).tip.yglob],...
                    [obj.panels(1).root.a,obj.panels(2).root.a,obj.panels(2).tip.a],yvec );
                % Geometric twist
                eps_vet  = interp1( [obj.panels(1).root.yglob,obj.panels(2).root.yglob,obj.panels(2).tip.yglob],...
                    [obj.panels(1).root.eps,obj.panels(2).root.eps,obj.panels(2).tip.eps],yvec );
                % Aerodynamic twist
                azl_vet  = interp1( [obj.panels(1).root.yglob,obj.panels(2).root.yglob,obj.panels(2).tip.yglob],...
                    [obj.panels(1).root.alpha0l,obj.panels(2).root.alpha0l,obj.panels(2).tip.alpha0l], yvec );
                % cm@ac for profiles
                cm0_vet  = interp1( [obj.panels(1).root.yglob,obj.panels(2).root.yglob,obj.panels(2).tip.yglob],...
                    [obj.panels(1).root.cmac,obj.panels(2).root.cmac,obj.panels(2).tip.cmac], yvec );
                % profile's aerodynamic centers
                xac_vet  = interp1( [obj.panels(1).root.yglob,obj.panels(2).root.yglob,obj.panels(2).tip.yglob],...
                    [obj.panels(1).root.x_ac,obj.panels(2).root.x_ac,obj.panels(2).tip.x_ac], yvec );
                XleMAC   = obj.meanprofile.xglob;
                MAC      = obj.meanprofile.c;
                %Cla      = obj.prf3DClean.a;
            end
            n_alpha = length( alpha );
            % Basic Load
            cm_add  = zeros( n_alpha,1 ); 
            cm0_int = cm0_vet(:).*cvet(:).^2;
            cm0     = obj.trapezoidal_int( cm0_int,yvec );
            cm0     = cm0*2/( obj.Sw*MAC);
            % Homemade do..while for cm
            tol = 1e-4; x_ref = 0.25;%cm_a = 0;
            % Iteration 0
            x1 = XleMAC+x_ref*MAC - ( x_vet(:) - x_vet(1) + cvet(:).*xac_vet(:) ) ;
            for i_alpha = 1:n_alpha
                cm_add_trap       = 0.5*( alpha(i_alpha)+eps_vet(:)-azl_vet(:) ).*cla_vet(:).*cvet(:).*x1;
                cm_add( i_alpha ) = obj.trapezoidal_int( cm_add_trap,yvec );
                cm_add( i_alpha ) = cm_add( i_alpha ) *2/( obj.Sw*MAC);
            end
            as     = polyfit( alpha,cm0+cm_add,1); cm_a = as(1); % Cm-alpha curve slope
            x_ref  = x_ref - cm_a/Cla;
            % Iterations
            while abs( cm_a ) > tol
                x1 = XleMAC+x_ref*MAC - ( x_vet(:) - x_vet(1) + cvet(:).*xac_vet(:) ) ;
                for i_alpha = 1:n_alpha
                    cm_add_trap       = ( alpha(i_alpha)+eps_vet(:)-azl_vet(:) ).*cla_vet(:).*cvet(:).*x1;
                    cm_add( i_alpha ) = obj.trapezoidal_int( cm_add_trap,yvec );
                    cm_add( i_alpha ) = cm_add( i_alpha )*2/( obj.Sw*MAC);
                end
                as     = polyfit( alpha,cm0+cm_add,1); cm_a = as(1); % Cm-alpha curve slope
                x_ref  = x_ref - cm_a/Cla;  % Updates the position of the wing ac
            end
        end
        
        function int = trapezoidal_int (~,f_vet,x_vet)
        %trapezoidal_int: function that calculates the integral for a
        %function f(x) using teh trapezoidal integration
        %INPUT
        %   f_vet: vector of sampled points of integral function
        %   x_vet: vector of sampling points for f
        nsec = length( x_vet ); nfun =  length( f_vet );
        if nsec ~= nfun
            error('Number of samples not equal')
        end
        int = 0;
        for i = 2:nsec
            int = int + 0.5*( f_vet(i) + f_vet(i-1) )*( x_vet(i)-x_vet(i-1) );
        end

        end
        
        function [varG,varA,coords] = HLAssign(obj,varG,nM,n_i)
            %HLAssign: ricava per interpolazione i valori delle
            %caratteristiche geometriche e aerodinamiche dei profili che
            %delimitano i flaps, cui posizione sull'apertura è contenuta
            %nella colonna 9 di varG
            %   varG: vettore contenente le dimensioni geometriche degli
            %       ipersostentatori
            %   nM: numero di condizioni di volo
            %   varA: vettore contenente le caratteristiche aerodinamiche
            %       dei profili che delimitano i flaps
            
            %nM = length(varA(:,1))*0.5;
            if nargin < 4
                varA = NaN(2*nM,8);
                n_i = 8; % it is used for retrocompatibility, because the original code needed less aerodynamic inputs
            else
               varA = NaN(2*nM,10); 
            end
            ny = 0; j=1;
            coords = NaN(2,2);
            while ny<2
                [varGtmp,ypts] = obj.panels(j).panelInterp(...
                    varG(:,9));

                if ~isempty(ypts) && ny == 0
                    varA(1:length(ypts)*nM,:) = obj.panels(j).panelInterp(varG(:,9),n_i); % a
                    varG(1:length(ypts),1:8) = varGtmp;
                    coords(1:length(ypts),1) = obj.panels(j).root.xglob +  (varG(1:length(ypts),9) - obj.panels(j).root.yglob)*tan( obj.panels(j).sweep*pi/180); %coordinata x
                    coords(1:length(ypts),2) = obj.panels(j).root.zglob +  (varG(1:length(ypts),9) - obj.panels(j).root.yglob)*tan( obj.panels(j).dihedral*pi/180); %coordinata z
                elseif ~isempty(ypts) %&& ny == 1
                    varA(nM+1:2*nM,:) = obj.panels(j).panelInterp(varG(:,9),n_i); % a
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
        function CL = lift_eval(~,alpha,profile)
            alpha_n = length( alpha ); CL = nan( alpha_n,1 );
            cf_mat = [1 profile.alphastar profile.alphastar^2 profile.alphastar^3;...
                1 profile.alphamax profile.alphamax^2 profile.alphamax^3;...
                0 1 2*profile.alphastar 3*profile.alphastar^2;...
                0 1 2*profile.alphamax 3*profile.alphamax^2];
            RHS   = [profile.clstar;profile.clmax;profile.a;0];
            cfs = cf_mat\RHS; % Coefficients for polinomial approximation of CL-alpha curve
            for i = 1:alpha_n
                if alpha(i) < min( profile.alpha0l,-1 )
                    % alpha< alpha_0L
                    % Fixes CL at CL alpha0L-1
                    CL(i) = profile.cl0 + profile.a*( profile.alpha0l-1 ); 
                elseif alpha(i) < profile.alphastar
                    % alpha0l < alpha < alpha*
                    % Linear section of CL
                    CL(i) = profile.cl0 + profile.a*alpha(i);
                elseif alpha(i) < profile.alphamax
                    % alpha* < alpha < alpha_max
                    % CL is given as a cubic polinomial
                    CL(i) = 0;
                    for k = 1:4
                        CL(i) = CL(i) + cfs(k)*alpha(i)^(k-1);
                    end
                else
                    % alpha > alpha max
                    % Fixes CL as CL_max
                    CL(i) = profile.clmax;
                end
            end
        end

        function profClass = aero3Dwing(obj, flg, M, deltaF, deltaS,Fcalc)
            %aero3Dwing calcola le caratteristiche aerodinamiche dell'ala
            %tridimensionale
            %INPUT
            %   profClass: oggetto classe profilo sul quale salvare i dati
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
                if isnan( Midx(nMi) )
                    error("Nessun dato trovato per il mach imposto")
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
                    obj.csocAvg           = csavg;
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
                        SflapTot   = SflapTot + obj.flaps(i).S;
                    end
                    for i = 1:obj.nslats
                        SslatTot   = SslatTot + obj.slats(i).S;
                        cbarocSlat = cbarocSlat + ...
                            (obj.slats(i).root.cextoc + obj.slats(i).tip.cextoc)/(obj.nslats*2);
                    end
                    obj.Sflaps    = SflapTot; obj.Sslats = SslatTot;
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
                        profClass.M(idx>0)         = NaN; profClass.a(idx>0) = NaN;
                        profClass.cl0(idx>0)       = NaN; profClass.clstar(idx>0) = NaN; 
                        profClass.clmax(idx>0)     = NaN; profClass.alphamax(idx>0)= NaN; 
                        profClass.alphastar(idx>0) = NaN;  profClass.alpha0l(idx>0)= NaN; 
                        profClass.cmac(idx>0)      = NaN; 
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

                obProf.clstar(Midx(n))  = obj.meanprofile.clstar(Midx(n));
                obProf.alpha0l(Midx(n)) = obj.meanprofile.alpha0l(Midx(n));
                obProf.cl0(Midx(n))     = obProf.a(Midx(n)) * (-obProf.alpha0l(Midx(n)));

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
                obj.flaps(i).root = obj.flaps(i).root.HLauxVariables('dCl02D', dclo);
                obj.flaps(i).root = obj.flaps(i).root.HLauxVariables('cbaroc', cboc);
                obj.flaps(i).root = obj.flaps(i).root.HLauxVariables('alphaDeltaf', alphadf);
                obj.flaps(i).root = obj.flaps(i).root.HLauxVariables('dClmax2D', dClMaxHLFun(deltaF,...
                    obj.flaps(i),'root'));
                [dclo,cboc,alphadf] = dCl02DHLFun(deltaF,obj.flaps(i),'tip',Midx);
                obj.flaps(i).tip = obj.flaps(i).tip.HLauxVariables('dCl02D', dclo);
                obj.flaps(i).tip = obj.flaps(i).tip.HLauxVariables('cbaroc', cboc);
                obj.flaps(i).tip = obj.flaps(i).tip.HLauxVariables('alphaDeltaf', alphadf);
                obj.flaps(i).tip = obj.flaps(i).tip.HLauxVariables('dClmax2D', dClMaxHLFun(deltaF,...
                    obj.flaps(i),'tip'));

                % Calcolo Pesi Corde Flaps
                k1 = 0.5*obj.flaps(i).root.c*obj.flaps(i).b/obj.flaps(i).S;
                k2 = 0.5*obj.flaps(i).tip.c*obj.flaps(i).b/obj.flaps(i).S;

                % Calcoli effetto flaps su cm
                weiS = k1 + k2 + weiS; % Somma dei pesi
                cfavg = cfavg + k1*obj.flaps(i).root.cfoc + k2*obj.flaps(i).tip.cfoc;
                %dClMax2D_flap
                %dClMax2D medio
                obj.flaps(i) = obj.flaps(i).HLauxVariables('dCLmax_mean', k1*obj.flaps(i).root.HLauxVariables('dClmax2D')+...
                    k2*obj.flaps(i).tip.HLauxVariables('dClmax2D'));

                sweep025 = ...
                    obj.sweepChange(obj.panels(1).sweep,0,0.25,obj.AR,obj.TR);
                Ksweep = (1-0.08*cos(sweep025*pi/180).^2)*cos(sweep025*pi/180).^(3/4) ;
                SfoS = 2*obj.flaps(i).S/obj.Sw;
                %dCLmax
                obj.flaps(i) = obj.flaps(i).coeffsShiftFun('dCLmax',SfoS*Ksweep*...
                    obj.flaps(i).HLauxVariables('dCLmax_mean'));

                %%dCl02D_flap
                %dCl02D medio
                obj.flaps(i) = obj.flaps(i).HLauxVariables('dCl02D_mean', ...
                    k1*obj.flaps(i).root.HLauxVariables('dCl02D')+...
                    k2*obj.flaps(i).tip.HLauxVariables('dCl02D'));
                % a medio flap
                obj.flaps(i) = obj.flaps(i).HLauxVariables('a_mean',k1*obj.flaps(i).root.a(Midx)+...
                    k2*obj.flaps(i).tip.a(Midx));
                % alphaDeltaf medio
                obj.flaps(i) = obj.flaps(i).HLauxVariables('alphaDeltaf' , ...
                    k1*obj.flaps(i).root.HLauxVariables('alphaDeltaf')+...
                    k2*obj.flaps(i).tip.HLauxVariables('alphaDeltaf'));
                %Kb flap
               obj.flaps(i)  = obj.flaps(i).HLauxVariables('Kb', kbFun(...
                    [2*obj.flaps(i).root.yglob/obj.bw,...
                    2*obj.flaps(i).tip.yglob/obj.bw],obj.TR));
                %Kc flap
                obj.flaps(i)  = obj.flaps(i).HLauxVariables('Kc', kcFun( ...
                    obj.flaps(i).HLauxVariables('alphaDeltaf')...
                    ,obj.AR,1));
                %dCL0 flaps
                obj.flaps(i)  = obj.flaps(i).coeffsShiftFun('dCL0',...
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
                obj.slats(i).root = obj.slats(i).root.HLauxVariables('etaMax',...
                    etaMaxFun(obj.slats(i).root.LERc/obj.slats(i).root.tc,1)); % DatoLER/t potrebbe essere sbagliato
                obj.slats(i).root = obj.slats(i).root.HLauxVariables('etaDelta',...
                    etaDeltaFun(deltaS,1));
                obj.slats(i).root = obj.slats(i).root.HLauxVariables('CloDs',...
                    clOdsSlat(obj.slats(i).root.cfoc));
                obj.slats(i).root = obj.slats(i).root.HLauxVariables('dClmax_slat',...
                    obj.slats(i).root.HLauxVariables('etaMax')*...
                    obj.slats(i).root.HLauxVariables('etaDelta')*...
                    obj.slats(i).root.HLauxVariables('CloDs')*...
                    deltaS*obj.slats(i).root.cextoc);
                % Sezione Slat di Estremità
                obj.slats(i).tip = obj.slats(i).tip.HLauxVariables('etaMax',...
                    etaMaxFun(obj.slats(i).tip.LERc/obj.slats(i).tip.tc,1));
                obj.slats(i).tip = obj.slats(i).tip.HLauxVariables('etaDelta',...
                    etaDeltaFun(deltaS,1));
                obj.slats(i).tip = obj.slats(i).tip.HLauxVariables('CloDs',...
                    clOdsSlat(obj.slats(i).tip.cfoc));
                obj.slats(i).tip = obj.slats(i).tip.HLauxVariables('dClmax_slat',...
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

                obj.slats(i) = obj.slats(i).HLauxVariables('dCmax2D_slat',...
                    k1*obj.slats(i).root.HLauxVariables('dClmax_slat') + ...
                    k2*obj.slats(i).tip.HLauxVariables('dClmax_slat') );
                sweep025 = ...
                    obj.sweepChange(obj.panels(1).sweep,0,0.25,obj.AR,obj.TR);
                Ksweep = (1-0.08*cos(sweep025*pi/180).^2)*cos(sweep025*pi/180).^(3/4) ;
                obj.slats(i) = obj.slats(i).coeffsShiftFun('dCLmax_slat',...
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
