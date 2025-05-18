classdef AirData_class
    
    properties 
        % - Nome
        Manufacturer
        Family
        Name
        Code
        Mark % marker per la grafica
        % - Crew
        npax
        npil
        ncrew
        % - Masse
        MRamp   %Ramp Mass
        MTOM    %Max.Take-Off Mass
        MLand   % Max Land Mass
        MZF     %Max. Zero-Fuel: ovvero il peso senza carburante ma con il max payload imbarcato
        MPM     %Max Payload Mass: limite a quanto payload si può imbarcare, sarebbe il payload delp Pto A
        MFW     %Max Fuel Weight: limite alla capacità dei serbatoi
        PatPBM  %Payload Mass at Point B: massa del payload quando si riempiono completamente i serbatoi (Pto B)
        DPM     %Design Payload Mass: massa del payload imbarcata nel pto di design (Pto D)
        DFW     %Design Fuel Weight: massa di carburante da imbarcare nel pto di design (Pto D)
        OEM     %Operating Empty Mass
        TFM     % Trapped Fuel Mass
        BEM     %Basic Empty Mass
        OIM     %Operating Items Mass
        EM      %Empty Mass: OEI - Wfuel - Wcrew quindi si suppone ci siano degli standard items medi
        % - Carburante
        fuel_cap %Capacità globale dei rerbatoi [L]
        un_fuel_cap % QUantità di fuel inutilizzabiile [L]
        nTank %numero di serbatoi
        % - Dimensioni Fusoliera
        fus_Length %Lunghezza [m]
        fus_Height %Altezza [m]
        fus_Width %Larghezza [m]
        fus_fitRatio
        fus_vol_paxcomp %Volume Fusoliera per aera passeggeri [m^3]
        fus_vol_cockpit %Volume del cockpit [m^3]
        fus_vol_unequipped %Volume di fusoliera non equipaggiata (senza sedili,etc) [m^3]
        fus_cc_vol % Volume utile dei cargo compartments [m^3]
        fus_nCC % Numero di cargo compartments
        
        wing
        horizontal
        vertical
        % - Undercarriage
        wheel_track
        wheelbase
        wheel_turn_rad
        main_wheel_D
        main_wheel_Width
        n_wheels
        % - Nacelles
        nacelle_length
        nacelle_width
        nacelle_height
        nacelle_num
        nacelles_xpos
        nacelles_ypos
        % - PRESTAZIONI
        
        % - T/O
        TO_ISA_ST_SL
        TO_ISA_15_SL
        TO_ISA_ST_5000
        TO_ISA_15_5000
        V2 % vel. sup. ostacolo (IAS) [kts]
        % - LND
        LND_ISA_ST_SL
        LND_ISA_15_SL
        LND_ISA_ST_5000
        LND_ISA_15_5000
        LND_V_app
        
        % - Crociera
        V_cr % velocità crociera TAS [m/s]
        M_cr % Mach crociera
        h_cr % quota crociera [m]
        e_cr % fattore di oswald
        Cd0_cr
        rho_cr % Densità alla quota di crociera modello ISA [kg/m^3]
        throttle_cr % Manetta in crociera
        W_cr % Peso a inizio crociera [N]
        % - Ranges
        Rmax_pay % Range max payload [km]
        Rdes  % Design range [km]
        Rmax_fuel % Range max fuel with payload [km]
        Rferry  % Ferry Range [km]
        % - Thrust
        T0 % max static thrust [N]
        Tmax_cont % max continous thrust [N]
        TSFC
        % - ROCs
        ROC % struttura contenente ROC_max [m/s] e quota [m] alla quale si ha il ROC
        % Derived Data
        ToW
        WoS
    end
    methods
        function obj = AirData_class(dataFileName)
            %% Lettura File Input
            f_id = fopen(dataFileName,'r');
            % Grafica di Lettura
            disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
            % TODO _ skippare righe 
            try
            % Lattura nome del Velivolo
            temp = fgetl(f_id); obj.Manufacturer = strtok(temp);
            temp = fgetl(f_id); obj.Family = strtok(temp);
            temp = fgetl(f_id); obj.Name = strtok(temp);
            temp = fgetl(f_id); obj.Code = strtok(temp);
            disp(['% NOW READING: ',dataFileName, '%'] )
            disp(['% Aircraft: ',obj.Family,' - ',obj.Name] )
            disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%') 
            
            % - Crew
            temp = fgetl(f_id);
            tag_nac = 'Crew	';
            if ~strcmp(temp,tag_nac)
                error('Expected Crew fields');
            end
            obj.npax = fscanf(f_id,'%f '); temp = fgetl(f_id);
            obj.npil = fscanf(f_id,'%f '); temp = fgetl(f_id);
            obj.ncrew = fscanf(f_id,'%f '); temp = fgetl(f_id);
            
            %- Masse
            temp = fgetl(f_id);
            tag_nac = 'Weights:	';
            if ~strcmp(temp,tag_nac)
                error('Expected Weights fields');
            end
            
            obj.MRamp = fscanf(f_id,'%f ');
            temp = fgetl(f_id);% read the rest of the line as dummy text
            obj.MTOM = fscanf(f_id,'%f '); temp = fgetl(f_id);
            obj.MLand = fscanf(f_id,'%f '); temp = fgetl(f_id);
            obj.MZF = fscanf(f_id,'%f '); temp = fgetl(f_id);
            obj.MPM = fscanf(f_id,'%f '); temp = fgetl(f_id);
            obj.PatPBM = fscanf(f_id,'%f '); temp = fgetl(f_id);
            obj.MFW = fscanf(f_id,'%f '); temp = fgetl(f_id);
            obj.DPM = fscanf(f_id,'%f '); temp = fgetl(f_id);
            obj.DFW = fscanf(f_id,'%f '); temp = fgetl(f_id);
            obj.OEM = fscanf(f_id,'%f '); temp = fgetl(f_id);
            obj.BEM = fscanf(f_id,'%f '); temp = fgetl(f_id);
            obj.EM = fscanf(f_id,'%f '); temp = fgetl(f_id);

            % TODO skippare righe
            
            %- Carburante
            temp = fgetl(f_id);
            tag_lit = 'Fuel_(litres)	';
            if ~strcmp(temp,tag_lit)
                error('Expected Fuel Tank Capacity fields');
            end
            temp = fgetl(f_id);
            tag_lit = 'Unusable(Lt)'; obj.fuel_cap = 0; obj.nTank = 0;
            while ~strcmp(temp,tag_lit)
                temp_fuel = obj.fuel_cap;
                obj.fuel_cap = obj.fuel_cap +  fscanf(f_id,'%f '); 
                obj.nTank = obj.nTank + 1;
                temp = fgetl(f_id);
            end
            obj.fuel_cap = temp_fuel; obj.nTank = obj.nTank - 1;
            
            if ~strcmp(temp,tag_lit)
                error('Expected Unusable Fuel Tank Capacity fields');
            end
            obj.un_fuel_cap = fscanf(f_id,'%f '); temp = fgetl(f_id);
            %- Dimensioni: Fusoliera
            temp = fgetl(f_id);
            tag_dim = 'DIMENSIONS	';
            if ~strcmp(temp,tag_dim)
                error('Expected Fuselage fields');
            end
            obj.fus_Length   = fscanf(f_id,'%f '); temp = fgetl(f_id);
            obj.fus_Height   = fscanf(f_id,'%f '); temp = fgetl(f_id);
            obj.fus_Width    = fscanf(f_id,'%f '); temp = fgetl(f_id);
            obj.fus_fitRatio = fscanf(f_id,'%f '); temp = fgetl(f_id);
            % Volumi di Fusoliera
            temp = fgetl(f_id);
            tag_fus2 = 'Fuselage Volumes  (m^3)	';
            if ~strcmp(temp,tag_fus2)
                error('Expected Fuselage Volumes fields');
            end
            obj.fus_vol_unequipped = fscanf(f_id,'%f '); temp = fgetl(f_id);
            obj.fus_vol_paxcomp    = fscanf(f_id,'%f '); temp = fgetl(f_id);
            obj.fus_vol_cockpit    = fscanf(f_id,'%f '); temp = fgetl(f_id);
            % Volumi Cargo Bay
            temp = fgetl(f_id);
            tag_fus3 = 'Cargo Compartment	';
            if ~strcmp(temp,tag_fus3)
                error('Expected Cargo Compartment fields');
            end
            tag_dim = 'Wing:	'; obj.fus_cc_vol = 0; obj.fus_nCC = 0;
            while ~strcmp(temp,tag_dim)
                temp_vol = obj.fus_cc_vol;
                obj.fus_cc_vol = obj.fus_cc_vol +  fscanf(f_id,'%f '); 
                obj.fus_nCC = obj.fus_nCC + 1;
                temp = fgetl(f_id);
            end
            obj.fus_cc_vol = temp_vol; obj.fus_nCC = obj.fus_nCC - 1;
            
            %- Dimensioni Ala
            tag_surf = 'Wing:	';
            if ~strcmp(temp,tag_surf)
                error('Expected Wing fields');
            end
            try
                obj.wing = obj.readVals(obj.fus_Width,f_id);
                obj.WoS = obj.MTOM/obj.wing.Sw;
            catch ERR
                warning( ['Error while calculating Horizontal geometry',newline,ERR.message]);
                obj.WoS = nan;
            end
            %- Dimensioni Piano Orizzontale
            temp = fgetl(f_id);
            tag_surf = 'Horizontal:	';
            if ~strcmp(temp,tag_surf)
                error('Expected Horizontal fields');
            end
            y2_rexp = fscanf(f_id,'%f '); temp = fgetl(f_id); 
            try
                obj.horizontal = obj.readVals(y2_rexp*2,f_id);
            catch ERR
                warning( ['Error while calculating Horizontal geometry',newline,ERR.message]);
            end
            % Temporary
            temp = fgetl(f_id);
            temp = fgetl(f_id);
            temp = fgetl(f_id);
            temp = fgetl(f_id);
            temp = fgetl(f_id);
            
            %- Dimensioni Piano Verticale
            temp = fgetl(f_id);
            tag_surf = 'Vertical:	';
            if ~strcmp(temp,tag_surf)
                error('Expected Vertical fields');
            end
            y2_rexp = fscanf(f_id,'%f '); temp = fgetl(f_id);
            try
                obj.vertical = obj.readVals(y2_rexp*2,f_id);
            catch ERR
                warning( ['Error while calculating vertical geometry',newline,ERR.message]);
            end
            % Temporary
            temp = fgetl(f_id);
            temp = fgetl(f_id);
            temp = fgetl(f_id);
            temp = fgetl(f_id);
            temp = fgetl(f_id);
            
            %- Dimensioni: Undercarriage
            temp = fgetl(f_id);
            tag_und = 'Undercarriage:	';
            if ~strcmp(temp,tag_und)
                error('Expected Undercarriage fields');
            end
            obj.wheel_track      = fscanf(f_id, '%f '); temp = fgetl(f_id);
            obj.wheelbase        = fscanf(f_id, '%f '); temp = fgetl(f_id);
            obj.wheel_turn_rad   = fscanf(f_id, '%f '); temp = fgetl(f_id);
            obj.main_wheel_D     = fscanf(f_id, '%f '); temp = fgetl(f_id);
            obj.main_wheel_Width = fscanf(f_id, '%f '); temp = fgetl(f_id);
            obj.n_wheels         = fscanf(f_id, '%f '); temp = fgetl(f_id);

            %- Dimensioni: Nacelles
            temp = fgetl(f_id);
            tag_nac = 'Nacelle:	';
            if ~strcmp(temp,tag_nac)
                error('Expected Nacelles fields');
            end
            obj.nacelle_length = fscanf(f_id, '%f '); temp = fgetl(f_id);
            obj.nacelle_width  = fscanf(f_id, '%f '); temp = fgetl(f_id);
            obj.nacelle_height = fscanf(f_id, '%f '); temp = fgetl(f_id);
            obj.nacelle_num    = fscanf(f_id, '%f '); temp = fgetl(f_id);
            obj.nacelles_xpos  = fscanf(f_id, '%f '); temp = fgetl(f_id);
            obj.nacelles_ypos  = fscanf(f_id, '%f '); temp = fgetl(f_id);
            
            %- PRESTAZIONI
            %- T/O
            temp = fgetl(f_id);
            tag_nac = 'Performances:	';
            if ~strcmp(temp,tag_nac)
                error('Expected Performances fields');
            end
            temp = fgetl(f_id);
            tag_nac = 'Take-Off	';
            if ~strcmp(temp,tag_nac)
                error('Expected Take-Off Performances fields');
            end
            obj.TO_ISA_ST_SL   = fscanf(f_id, '%f '); temp = fgetl(f_id);
            obj.TO_ISA_15_SL   = fscanf(f_id, '%f '); temp = fgetl(f_id);
            obj.TO_ISA_ST_5000 = fscanf(f_id, '%f '); temp = fgetl(f_id);
            obj.TO_ISA_15_5000 = fscanf(f_id, '%f '); temp = fgetl(f_id);
            obj.V2 = fscanf(f_id, '%f '); temp = fgetl(f_id); 
            %- LND
            temp = fgetl(f_id);
            tag_nac = 'Landing	';
            if ~strcmp(temp,tag_nac)
                error('Expected Landing Performances fields');
            end
            obj.LND_ISA_ST_SL = fscanf(f_id, '%f '); temp = fgetl(f_id);
            obj.LND_ISA_15_SL  = fscanf(f_id, '%f '); temp = fgetl(f_id);
            obj.LND_ISA_ST_5000 = fscanf(f_id, '%f '); temp = fgetl(f_id);
            obj.LND_ISA_15_5000    = fscanf(f_id, '%f '); temp = fgetl(f_id);
            obj.LND_V_app    = fscanf(f_id, '%f '); temp = fgetl(f_id);
            %- Cruise
            temp = fgetl(f_id);
            tag_nac = 'Cruise	';
            if ~strcmp(temp,tag_nac)
                error('Expected Cruise Performances fields');
            end
            temp   = fscanf(f_id, '%f '); obj.V_cr = convvel(temp,'kts','m/s'); temp = fgetl(f_id);
            obj.M_cr   = fscanf(f_id, '%f '); temp = fgetl(f_id);
            temp   = fscanf(f_id, '%f ');obj.h_cr = convlength(temp,'ft','m'); temp = fgetl(f_id);
            obj.e_cr   = fscanf(f_id, '%f '); temp = fgetl(f_id);
            obj.Cd0_cr = fscanf(f_id, '%f '); temp = fgetl(f_id);
            % - Ranges
            temp = fgetl(f_id);
            tag_nac = 'Range	';
            if ~strcmp(temp,tag_nac)
                error('Expected Ranges fields');
            end
            temp   = fscanf(f_id, '%f '); obj.Rmax_pay = convlength(temp,'naut mi','km'); temp = fgetl(f_id);
            temp   = fscanf(f_id, '%f '); obj.Rdes = convlength(temp,'naut mi','km'); temp = fgetl(f_id);
            temp   = fscanf(f_id, '%f '); obj.Rmax_fuel = convlength(temp,'naut mi','km'); temp = fgetl(f_id);
            temp   = fscanf(f_id, '%f '); obj.Rferry = convlength(temp,'naut mi','km'); temp = fgetl(f_id);

            % - ROCs
            temp = fgetl(f_id);
            tag_nac = 'ROCs	';
            if ~strcmp(temp,tag_nac)
                error('Expected Rate-of-Climb fields');
            end
            
            tag_dim = 'Thrust	'; %obj.fus_cc_vol = 0; obj.fus_nCC = 0;
            idx = 1; temp = fscanf(f_id,'%f ');
            while ~strcmp(temp,tag_dim) && ~isempty( temp )
                obj.ROC.h(idx) = temp*0.305;
                temp = fgetl(f_id);
                obj.ROC.roc(idx) = fscanf(f_id,'%f ')*(0.305/60);
                temp = fgetl(f_id);
                idx = idx + 1;
                temp = fscanf(f_id,'%f ');
            end
            
            % - Spinta
            temp = fgetl(f_id);
            tag_nac = 'Thrust	';
            if ~strcmp(temp,tag_nac)
                error('Expected Thrust fields');
            end
            obj.T0        = fscanf(f_id, '%f '); temp = fgetl(f_id);
            obj.Tmax_cont = fscanf(f_id, '%f '); temp = fgetl(f_id);
            catch ERR
                warning( ['Reading from file failed ',newline,ERR.message] )
            end
            fclose(f_id);
            
            %% Calcoli Derivati
            %_Pesi
            obj = obj.weightEst;
            
            %Estrazione Aerodinamica
            obj.W_cr = 0.992*0.996*0.996*0.990*obj.MTOM*9.81; % Peso inizio crociera
            [~, ~, ~, obj.rho_cr] = atmosisa(obj.h_cr); % Densità inizio atmosfera
            
            % Atro
            obj.Mark = obj.marker_assign();
            
            if isempty(obj.T0)
                obj.ToW = nan;
            else
                obj.ToW = (obj.nacelle_num*obj.T0/9.81)/obj.MTOM;
            end
        end
        function obj = weightEst(obj)
            %weightEst: funzione che ricava i pesi dai pesi definiti.
            if isempty( obj.MTOM )
                % TODO
            else
                if ~isempty( obj.fuel_cap )
                    %                     fuel_cap in [Litri]
                    rho_f = 0.785; %[Kg/L]
                    check_temp = obj.fuel_cap * rho_f;
                end
                if isempty( obj.MFW )
                    obj.MFW = check_temp;
                else
                    disp( ' Max Fuel Weight assigned vs Calculaed from Tanks Capacity ' );
                    obj.checkLog( obj.MFW,check_temp,'MFW' );
                end

                if isempty( obj.MPM )&& isempty( obj.OEM )
                    warining('Unable to calculate Operative Empty Mass');
                    obj.OEM = nan;
                else
                    if isempty( obj.MPM )
                        obj.MPM = obj.MZF - obj.OEM;
                    else
                        % DEBUG
                        disp('Check pax')
                        ( 1 - obj.npax*215/(2.2046 * obj.MPM) ) *100
                        
                        check_temp = obj.MZF - obj.MPM;
                        if isempty( obj.OEM )
                            obj.OEM = check_temp;
                        else
                            obj.checkLog( obj.OEM,check_temp,'OEM' );
                        end
                    end
                end
                % ???? Check
                if isempty( obj.npax ) || isempty( obj.npil )
                    if isempty( obj.EM )
                        warning(' Cannot estimate Empty Weight' );
                        disp(' The following results will be WRONG' );
                        obj.EM = 0;
                    end
                else
                    % Calcolo della Massa della Crew. Il numero di
                    % assistenti di volo si ottiene dalla legge 1 ogni 50
                    % passeggeri
                    Mcrew = obj.npil*94 ;
                    ncrew_est = ceil( obj.npax/50 );
                    if isempty( obj.ncrew )
                        obj.ncrew = ncrew_est;
                    else
                        obj.checkLog( obj.ncrew,ncrew_est,'Ncrew' );
                    end
                    Mcrew = Mcrew + obj.ncrew*94 ;
                    % Calcolo della Massa del Trapped Fuel
                    if obj.MTOM < 1e5*0.4536
                        Mtf = 0;
                    else
                        Mtf = 0.005*obj.MTOM;
                    end
                    check_temp = obj.OEM - Mcrew - Mtf;
                    if isempty( obj.EM )
                        obj.EM = check_temp;
                    else
                        obj.checkLog( obj.EM,check_temp,'Ncrew' );
                        if obj.EM > obj.OEM
                            error(' Basic Empty Mass assigned is greater than Operating Empty Mass, Check Data' );
                        end
                    end
                end
                %obj.OIM = obj.OEM - obj.BEM;
                %check_temp = obj.OIM/obj.MTOM;
                %if check_temp>0.03 && check_temp < 0.06
                    disp( ['ME/MTOM is : ',num2str(check_temp*100)] );
%                 else
%                     warning( 'Operating Weigths fraction is unusual' );
%                     check_temp = 0.045;
%                     if obj.OEM - obj.MTOM*check_temp > 0
%                         disp( [' Assuming OIM as ',num2str(check_temp*100),'% of MTOM'] );
%                         obj.OIM = obj.MTOM*check_temp;
%                         obj.BEM = obj.OEM - obj.OIM;
%                     end
                %end
            end
        end
        
        function [lsurf_str,S,AR,TR,sweepLE,sweep025,cmac] = prepGeom(obj,lsurf_str,b,S,AR,TR,sweepLE_in)
            %prepGeom: function that arranges data readed from txt file
            %into variables used to define WingClass
            % INPUT:
            %   lsufr_str: struct contenente come campi i parametri
            %       geometrici della superficie%    
            %   b: aperura alare
            %   S: superficie alare
            TRi = 0;
            if isempty( lsurf_str.ckink )
                TRi = 1; lsurf_str.y_kink = 0;
                chk = 1; % flag to check if the wing hasn't got a kink
            else
                if isempty( lsurf_str.y_kink )
                    warning( 'Kink Position Required, cannot proceed' );
                    lsurf_str.y_kink = nan;
                    chk = 0;
                end
            end
            %chk  = TRi == 1;
            if isempty( lsurf_str.croot ) && isempty( lsurf_str.ctip )
               warning( 'Too few data, cannot proceed to geometric calculations' );
               lsurf_str.croot = nan; lsurf_str.ctip = nan;
            else
                if isempty( lsurf_str.ctip )
                    if isempty( TR )
                        warning( 'Cannot calculate ctip' );
                        lsurf_str.ctip = nan;
                    else
                        lsurf_str.ctip = TR*lsurf_str.croot;
                    end
                else
                    if isempty( lsurf_str.croot_exp ) || isempty( lsurf_str.y_rootexp ) || isempty( b )
                        warning( 'Cannot calculate croot' );
                        lsurf_str.croot = nan;
                    else
%                         chk  = TRi == 1;
                        if chk 
                            % no kink
                            lsurf_str.x_apex = lsurf_str.xle_exp + (lsurf_str.xle_tip - lsurf_str.xle_exp)/(0.5*b - lsurf_str.y_rootexp)*(0-lsurf_str.y_rootexp);
                            xte_root = lsurf_str.xle_exp + lsurf_str.croot_exp + ...
                                ( lsurf_str.xle_tip+lsurf_str.ctip - (lsurf_str.xle_exp+lsurf_str.croot_exp) )/(0.5*b - lsurf_str.y_rootexp)*(0-lsurf_str.y_rootexp);
                        else
                            % kink
                            lsurf_str.x_apex = lsurf_str.xle_exp + (lsurf_str.xle_kink -lsurf_str.xle_exp)/(lsurf_str.y_kink - lsurf_str.y_rootexp)*(0-lsurf_str.y_rootexp);
                             xte_root = lsurf_str.xle_exp + lsurf_str.croot_exp + ...
                                 ( lsurf_str.xle_kink + lsurf_str.ckink - (lsurf_str.xle_exp+lsurf_str.croot_exp) )/(lsurf_str.y_kink - lsurf_str.y_rootexp)*(0-lsurf_str.y_rootexp);
                        end
                        %xte_root = lsurf_str.xle_exp + lsurf_str.croot_exp + ...
                            %( lsurf_str.xle_kink + lsurf_str.ckink - (lsurf_str.xle_exp+lsurf_str.croot_exp) )/(lsurf_str.y_kink - lsurf_str.y_rootexp)*(0-lsurf_str.y_rootexp)*(1 - chk) + ... %interp con kink
                            %( lsurf_str.xle_tip+lsurf_str.ctip - (lsurf_str.xle_exp+lsurf_str.croot_exp) )/(0.5*b - lsurf_str.y_rootexp)*(0-lsurf_str.y_rootexp)*(chk); % interp con tip
                        lsurf_str.croot = xte_root - lsurf_str.x_apex;
                    end
                end
                if chk
                    tmp = 0.5*b*(lsurf_str.croot+lsurf_str.ctip);
                else
                    tmp = lsurf_str.y_kink*(lsurf_str.croot+lsurf_str.ckink) ...
                        + (0.5*b - lsurf_str.y_kink)*(lsurf_str.ckink+lsurf_str.ctip);
                end
                if isempty( S )
                    S = tmp;
                else
                    obj.checkLog( S,tmp,'S' );
                end
                tmp = b^2/S;
                if isempty( AR )
                    AR = tmp;
                else
                    obj.checkLog( AR,tmp,'AR' );
                end
                tmp = lsurf_str.ctip/lsurf_str.croot;
                if isempty( TR )
                    TR = tmp;
                else
                    obj.checkLog( TR,tmp,'TR' );
                end
                
                [sweepLE,~] = obj.eqwing_fun( sweepLE_in,lsurf_str,S,chk ); 
                lsurf_str.toc_avg = obj.weightendavg( lsurf_str,S,chk  );
                [cmac,ymac] = obj.macalc( lsurf_str,TR,chk);
                sweep025 = obj.sweepChange(S,0,0.25,AR,TR);
            end
        end
        
        function ls_obj = readVals(obj,d_fus,f_id)
            %readVals: function that reads from a txt file the geometric
            %parametres for a lifting surface
            temp_sects = nan(3,4); %geom_sects_aus = nan(2,1);
            for i = 1:2
                % i = 1
                % Reads c_root, c_root_exp c_kink c_tip
                % i = 2
                % Reads x_apex, xLE_root_exp xLE_kink xLE_tip
                for j = 1:4
                    temp = fscanf(f_id, '%f ');
                    if isempty( temp )
                        temp = nan;
                    end
                    temp_sects(i,j) = temp;
                    temp = fgetl(f_id);
                end
            end
            temp_sects      = temp_sects';
            if isempty( d_fus )
                % if the exposed root is not defined y_kink is set to nan
                temp_sects(2,3) = nan;
            else
                temp_sects(2,3) = d_fus;            % 2*y_cexp
            end
            for i = 1:2
                % Reads y_kink y_tip
                temp = fscanf(f_id, '%f '); 
                if isempty( temp )
                    temp = nan;
                end
                temp_sects(i+2,3) = temp;
                temp = fgetl(f_id);
            end
            temp_sects(3,3) = 2*temp_sects(3,3); % 2*y_kink
            %geom_sects_aus = geom_sects_aus(:)';
            %S        = fscanf(f_id, '%f '); 
            temp = fgetl(f_id); % skips S read
            
            temp = fgetl(f_id);     % skips sweep read sweepLE  = fscanf(f_id, '%f '); 
            temp = fgetl(f_id);     % skips sweep@0.25 read sweep025 = fscanf(f_id, '%f '); 
            dihedral = fscanf(f_id, '%f '); temp = fgetl(f_id);
            
            temp = fgetl(f_id);     % skips AR read AR = fscanf(f_id, '%f ');
            temp = fgetl(f_id);     % skips TR read TR = fscanf(f_id, '%f ');
            % Max thickness distribution
            toc_vett = nan(4,1);
            for i=1:4
                % Reads toc root,kink,tip and avg
                temp = fscanf(f_id, '%f ');
                if isempty( temp)
                    temp = nan;
                end
                toc_vett(i) =  temp;
                temp = fgetl(f_id);
            end

            i_ang = fscanf(f_id, '%f ');
            if isempty( i_ang )
                i_ang = 0;
            end
            temp = fgetl(f_id);
            
            [bs,sweeps,apexC,sectsGeom,dihedrals] = obj.checkGeom(temp_sects,toc_vett,dihedral);
            M = nan; sectsAero = nan( length(bs)+1,8 );
            % Definition of Lifting Surface Class Object
            ls_obj = WingClass(bs,sweeps,dihedrals,i_ang,apexC,M,...
                sectsGeom, sectsAero);
        end

        function [bvet,sweep_vect,apexC,sect_vec,dih_vec] = checkGeom(obj,temp_sects,toc_vec,dih_vec)
            %checkGeom: function that checks if there are enough
            %geometrical data to define the wing and prepares the input
            %variables for the wingClass constructor
            %INPUT
            %   temp_sects: 4X3 matrix containing chords, x_LE coords,
            %       2*y_sect for root, exp. root, kink and tip sections
            %   toc_vec: vector containing the t/c_max at each section plus
            %       the averange value.
            %   dih_vec: vector with dihedral angles [deg] for each panel
            if ( isnan( temp_sects(1,1) )&& isnan( temp_sects(1,2) ) ||  isnan( temp_sects(2,1) )&& isnan( temp_sects(2,2) ) ) ... % isnan c_root && x_apex || isnan croot_exp %% df/2
                    && isnan( temp_sects(4,1) ) && isnan( temp_sects(4,2) ) % isnan(ctip) && isnan(b)
                warning('Not enough data for computing wing geometry');
            else
                if xor( isnan(temp_sects(3,1)),isnan(temp_sects(3,2)) ) % XOR( c_kink,2*y_kink )
                    warning('Incomplete kink data, cannot proceed');
                else
                    apexC = zeros(3,1);             % Init. wing apex vector
                    
                    if isnan( temp_sects(3,1) )     % isnan ckink
                        sect_vec = nan(2,8);        % Init. the geometric params. for the profile class
                        % single panel wing
                        if isnan( temp_sects(1,1) ) % isnan croot
                            % defined using cexp_root
                            [temp1,temp2,sweep_vect] = obj.extrapolate_root( temp_sects(2,1),temp_sects(2,2),0.5*temp_sects(2,3),...
                                temp_sects(4,1),temp_sects(4,2),0.5*temp_sects(4,3) );
                            
                            sect_vec(1,1) = temp1;  % Assign croot
                            apexC(1)      = temp2;  % Assign wing apex (croot xLE)
                            
                        else
                            % defined using croot
                            sweep_vect = atan( (temp_sects(4,2)-temp_sects(1,2)) / (0.5*temp_sects(4,3)) )*180/pi; % atan( (x_tip-x_apex) / b/2 )
                            
                            sect_vec(1,1) = temp_sects(1,1);
                            apexC(1)      = temp_sects(1,2);
                        end
                        sect_vec(2,1) = temp_sects(4,1);
                        bvet          = 0.5*temp_sects(4,3); % Assign panel span
                        % Check on input dihedrals
                        if length(dih_vec) > 1 || isempty(dih_vec)
                            warning( 'There are too many dihedrals for the wing panels considered' );
                            dih_vec = nan;
                        end
                        j = 1;
                        for i = [1,3]
                            if isnan( toc_vec(i) )
                                if ~isnan( toc_vec(4) )
                                    % assign t/c_avg to every section
                                    sect_vec(1:2,3) = toc_vec(4)*ones(2,1);
                                end
                                % no informations about t/c
                                break
                            end
                            sect_vec(j,3) = toc_vec(i);
                            j = j+1;
                        end
                        
                    else
                        %double paneled wing
                        sect_vec      = nan(3,8); % Init. the geometric params. for the profile class
                        sweep_vect    = nan(2,1);
                        sweep_vect(2) = atan( 2*(temp_sects(4,2)-temp_sects(3,2)) / ( temp_sects(4,3)-temp_sects(3,3) ) )*180/pi; %Sweep_out = atan( (x_tip-x_kink)/( b/2 - y_kink) )
                        bvet          = nan(2,1); % Init. panel span vector
                        if isnan( temp_sects(1,1) ) % isnan croot
                            % defined using cexp_root
                            [temp1,temp2,temp3] = obj.extrapolate_root( temp_sects(2,1),temp_sects(2,2),0.5*temp_sects(2,3),...
                                temp_sects(3,1),temp_sects(3,2),0.5*temp_sects(3,3) );
                            sweep_vect(1) = temp3;
                            apexC(1)      = temp2;
                            sect_vec(1,1) = temp1;
                        else
                            % defined using c_root
                            sweep_vect(1) = 2*(temp_sects(3,2)-temp_sects(1,2)) / ( temp_sects(3,3)-temp_sects(1,3)); %atan( (x_kink-x_apex)/( y_kink) )
                            apexC(1)      = temp_sects(1,2);                % x_apex
                            sect_vec(1,1) = temp_sects(1,1);                % croot
                            sweep_vect(1) = atan( sweep_vect(1) )*180/pi;   % Sweep in [deg]
                        end
                        sect_vec(2:3,1) = temp_sects(3:4,1);                        % ckink ctip
                        bvet(1,1)       = 0.5*temp_sects(3,3);                      % assign inner panel span
                        bvet(2,1)       = 0.5*( temp_sects(4,3) - temp_sects(3,3) );% assign outer panel span
                        % Check on Dihedrals
                        
                        switch length( dih_vec )
                            case 0
                                warning( 'Dihedral not defined' );
                                dih_vec = nan(2,1);
                            case 1
                                % Assumes same dihedral for both panels
                                dih_vec = [1,1]*dih_vec;
                            case 2
                                dih_vec = dih_vec(:)';
                            otherwise
                                warning( 'There are too many dihedrals for the wing panels considered' );
                        end
                        
                        % Sweep Continuity check
                        if abs( 1-sweep_vect(1)/sweep_vect(2) ) < 5e-3
                            sweep_vect(1) = ( sweep_vect(1) +sweep_vect(2) )*0.5;
                            sweep_vect(2) = sweep_vect(1);
                        end
                        j = 1;
                        for i = 1:3
                            if isnan( toc_vec(i) )
                                if ~isnan( toc_vec(4) )
                                    % assign t/c_avg to every section
                                    sect_vec(1:3,3) = toc_vec(4)*ones(2,1);
                                end
                                % no informations about t/c
                                break
                            end
                            sect_vec(j,3) = toc_vec(i);
                            j = j+1;
                        end
                        
                    end
                end
            end
            
        end
        
        function [croot,x_apex,tan_sweep] = extrapolate_root(~,c1,x1,y1,c2,x2,y2)
            %extrapolate_root: function that given two sections
            %extrapolates the coordinates and chord of the root
            %INPUT
            %   c1,x1,y1: chord and coordinates of the most inboard
            %       section
            %   c2,x2,y2: chord and coordinates of the most outboard
            %       section
            %OUTPUTS:
            %   tan_sweep: sweep angle in [deg]
            tan_sweep = (x2-x1) / ( y2-y1 ) ; % (x_tip-xLE_exp) / (b/2 - df/2) )
            x_apex    = x1 - tan_sweep*y1; %x_apex = xLE_croot_ecp - tan(sweep_in)*df/2
            tan_TE    = (x2+c2-x1-c1) / ( y2-y1 );
            temp      = x1+c1 - tan_TE*y1; %x_apex = xTE_croot_ecp - tan(sweep_in)*df/2
            croot = temp - x_apex;
            tan_sweep = atan( tan_sweep )*180/pi;
        end
        
        function [sweepLE_equiv,croot_equiv] = eqsurf( ~,lsurf_str,S,kinkFLG )
            
            sweepLE_equiv = 2*lsurf_str.y_tip*(lsurf_str.xle_tip - lsurf_str.x_apex);
            croot_equiv = S/lsurf_str.y_tip - lsurf_str.ctip;
            if kinkFLG
                % no kink
                sweepLE_equiv = atan( sweepLE_equiv/( 0.5*(2*lsurf_str.y_tip)^2 ) );
            else
                % kink
                tmp = ( lsurf_str.x_apex + lsurf_str.xle_kink - 2*lsurf_str.x_apex )*lsurf_str.y_kink;
                tmp = tmp + ( lsurf_str.xle_kink + lsurf_str.xle_tip - 2*lsurf_str.x_apex )*(lsurf_str.y_tip - lsurf_str.y_kink);
                sweepLE_equiv = atan( ( sweepLE_equiv-tmp )/( 0.5*(2*lsurf_str.y_tip)^2 ) );%( lsurf_str.y_tip^2 ) );
            end
           sweepLE_equiv = sweepLE_equiv*180/pi;
            
        end
        
        function checkLog( ~,input2check,calculated_value,nam,tol )
           %checkLog: funzione che scrive in un file il confronto tra i due valori in input
           if nargin < 5
               tol = 1e-1;
           end
           err = abs( 1 - calculated_value/input2check );
           if err > tol
               warning( [nam,' calculated has an error that excedes ',num2str(tol),'% from ref. value']);
           end
            
        end
        
        function avg = weightendavg( ~,lsurf_str,S,kinkFLG )
            %weightendavg: funzione che calcola le medie pesate cone le
            %aree delle sezioni delle ali
            Kroot = lsurf_str.croot*lsurf_str.y_kink/S; Ktip = lsurf_str.ctip*(lsurf_str.y_tip - lsurf_str.y_kink)/S;
            if kinkFLG
                Kkink = ( lsurf_str.croot*lsurf_str.y_kink+lsurf_str.ckink*(lsurf_str.y_tip - lsurf_str.y_kink) )/S;
                avg = lsurf_str.toc_kink*Kroot + ...
                    lsurf_str.toc_kink*Kkink + lsurf_str.toc_tip*Ktip;
            else
                 avg = lsurf_str.toc_kink*Kroot + lsurf_str.toc_tip*Ktip;
            end
        end
        
        function [sweepLE_temp,croot_eq] = eqwing_fun( obj,sweepLE,lsurf_str,S,kinkFLG )
            [sweepLE_temp,croot_eq] = obj.eqsurf( lsurf_str,S,kinkFLG );
            if ~isempty( sweepLE )
                obj.checkLog( sweepLE,sweepLE_temp,'Sweep at LE' );
            else
                
            end
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
        
        function [cmac,ymac] = macalc( ~,lsurf_str,TR,chk )
            if chk
                % Ala senza Kink
                cmac = (2/3)*lsurf_str.croot*(1+TR+TR^2)/(1+TR);
                ymac = lsurf_str.y_tip/3 * (1+2*TR)/(1+TR);
            else
                % Ala con Kink
                TRi = lsurf_str.ckink/lsurf_str.croot;
                A1 = lsurf_str.y_kink / (lsurf_str.y_tip - lsurf_str.y_kink);
                A2 = ( TRi+TR + A1 * (1+TRi) );
                cmac = (2/3)*lsurf_str.croot* ( (TRi^2 + TRi*TR + TR^2) ...
                    + A1*( 1 + TRi + TRi^2 ) );
                cmac = cmac/A2;
                ymac = (lsurf_str.y_tip - lsurf_str.y_kink)/3  * TRi*(1+2*TR*TRi) + ...
                    A1*(2+TRi*TR) + A1^2*(1+2*TRi);
                ymac = ymac/A2;
            end
        end
        
        function RES = check_delimiter(inp,namecheck)
            if ~strcmp(inp,namecheck)
                error('Expected Fuselage fields');
            end
        end
        
        function mark = marker_assign( obj )
        %marker_assign: assegna un marker in base alla famiglia di aereo
        switch obj.Manufacturer
            case 'Boeing'
                mark = 'o';
            case 'Airbus'
                mark = 'd';
            case 'Embraer'
                mark = 's';
            case 'Comac'
                mark = 'p';
            otherwise
                mark = '^';
        end
            
        end
    
    end
    
    
    
end
