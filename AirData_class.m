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
        % - Dimensioni Ala
        wing %struttura che ocntiene come campi wing.c = corde delle sezioni wing.x,wing.y = coordinate x/y delle sezioni wing.toc = t/c delle sezioni + t/c medio
        % 1 = root 2 = root_exp 3 = kink 4 = tip
        wing_b %apertura alare [m]
        wing_S %superficie alare
        wing_sweepLE %sweep al Leading edge [deg]
        wing_025sweep %sweep al 25 della corda
        wing_dihedral %Diedro [deg]
        wing_AR %Aspect Ratio
        wing_TR %Taper Ratio
        wing_toc_avg %t/c medio
        wing_mac %corda media aerodinamica [m]
        % - Dimensioni Piano Orizzontale
        hor %struttura che ocntiene come campi wing.c = corde delle sezioni wing.x,wing.y = coordinate x/y delle sezioni
        % 1 = root 2 = root_exp 3 = kink 4 = tip
        hor_b %apertura alare [m]
        hor_S %superficie alare [m^2]
        hor_sweepLE
        hor_sweep025
        hor_dihedral
        hor_AR 
        hor_TR %taper ratio
        hor_mac %corda media aerodinamica [m]
        hor_L %distanza tra i centri aerodinamici di ala e piano [m]
        hor_ShoS %rapporto S_hor/S_wing
        hor_VolRat %rapporto volumetrico piano di coda
        hor_cel %corda dell'elevator
        hor_ceoc % c_el/mac_hor
        % - Dimensioni Piano Verticale
        vert
        vert_b
        vert_S
        vert_sweepLE
        vert_sweep025
        vert_dihedral
        vert_AR
        vert_TR
        vert_mac
        vert_L
        vert_ShoS
        vert_VolRat
        vert_crud %Corda del timone
        vert_croc %c_rud/mac_ver
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
            tag_wing = 'Wing:	';
            if ~strcmp(temp,tag_wing)
                error('Expected Wing fields');
            end
            obj.wing.croot    = fscanf(f_id, '%f '); temp = fgetl(f_id); %croot
            obj.wing.croot_exp= fscanf(f_id, '%f '); temp = fgetl(f_id); %croot_exp
            obj.wing.ckink    = fscanf(f_id, '%f '); temp = fgetl(f_id); %ckink
            obj.wing.ctip     = fscanf(f_id, '%f '); temp = fgetl(f_id); %ctip
            obj.wing.x_apex   = nan; % wing apex
            obj.wing.xle_exp  = fscanf(f_id, '%f '); temp = fgetl(f_id); %xle_cexp
            obj.wing.xle_kink = fscanf(f_id, '%f '); temp = fgetl(f_id); %xle_kink
            obj.wing.xle_tip  = fscanf(f_id, '%f '); temp = fgetl(f_id); %xle_tip
            obj.wing.y_root = 0; obj.wing.y_rootexp = obj.fus_Width*0.5;
            obj.wing.y_kink   = fscanf(f_id, '%f '); temp = fgetl(f_id); %y_link
            obj.wing_b        = fscanf(f_id, '%f '); temp = fgetl(f_id);
            obj.wing.y_tip = obj.wing_b*0.5;
            obj.wing_S        = fscanf(f_id, '%f '); temp = fgetl(f_id);
            obj.wing_sweepLE  = fscanf(f_id, '%f '); temp = fgetl(f_id);
            obj.wing_025sweep = fscanf(f_id, '%f '); temp = fgetl(f_id);
            obj.wing_dihedral = fscanf(f_id, '%f '); temp = fgetl(f_id);
            obj.wing_AR       = fscanf(f_id, '%f '); temp = fgetl(f_id);
            obj.wing_TR       = fscanf(f_id, '%f '); temp = fgetl(f_id);
            obj.wing.toc_root = fscanf(f_id, '%f '); temp = fgetl(f_id);
            obj.wing.toc_kink = fscanf(f_id, '%f '); temp = fgetl(f_id);
            obj.wing.toc_tip  = fscanf(f_id, '%f '); temp = fgetl(f_id);
            obj.wing.toc_avg  = fscanf(f_id, '%f '); temp = fgetl(f_id); %sarebbe t/c_avg
            obj.wing_mac      = fscanf(f_id, '%f '); temp = fgetl(f_id);

            %- Dimensioni Piano Orizzontale
            temp = fgetl(f_id);
            tag_hor = 'Horizontal:	';
            if ~strcmp(temp,tag_hor)
                error('Expected Horizontal fields');
            end
            
            obj.hor.croot     = fscanf(f_id, '%f '); temp = fgetl(f_id); %croot
            obj.hor.croot_exp = fscanf(f_id, '%f '); temp = fgetl(f_id); %croot_exp
            obj.hor.ckink     = fscanf(f_id, '%f '); temp = fgetl(f_id); %ckink
            obj.hor.ctip      = fscanf(f_id, '%f '); temp = fgetl(f_id); %ctip
            obj.hor.x_apex    = nan;
            obj.hor.xle_exp   = fscanf(f_id, '%f '); temp = fgetl(f_id); %xle_exp
            obj.hor.xle_kink  = fscanf(f_id, '%f '); temp = fgetl(f_id); %xle_kink
            obj.hor.xle_tip   = fscanf(f_id, '%f '); temp = fgetl(f_id); %xle_tip
            obj.hor.y_root    = 0;
            obj.hor.y_rootexp = fscanf(f_id, '%f '); temp = fgetl(f_id); %y_root_exp
            obj.hor.y_kink    = fscanf(f_id, '%f '); temp = fgetl(f_id); %y_kink
            obj.hor_b         = fscanf(f_id, '%f '); temp = fgetl(f_id);
            obj.hor.y_tip     = obj.hor_b * 0.5;
            obj.hor_S         = fscanf(f_id, '%f '); temp = fgetl(f_id);
            obj.hor_sweepLE   = fscanf(f_id, '%f '); temp = fgetl(f_id);
            obj.hor_sweep025  = fscanf(f_id, '%f '); temp = fgetl(f_id);
            obj.hor_dihedral  = fscanf(f_id, '%f '); temp = fgetl(f_id);
            obj.hor_AR        = fscanf(f_id, '%f '); temp = fgetl(f_id);
            obj.hor_TR        = fscanf(f_id, '%f '); temp = fgetl(f_id);
            obj.hor.toc_root  = fscanf(f_id, '%f '); temp = fgetl(f_id);
            obj.hor.toc_kink  = fscanf(f_id, '%f '); temp = fgetl(f_id);
            obj.hor.toc_tip   = fscanf(f_id, '%f '); temp = fgetl(f_id);
            obj.hor.toc_avg   = fscanf(f_id, '%f '); temp = fgetl(f_id); %sarebbe t/c_avg
            obj.hor_mac       = fscanf(f_id, '%f '); temp = fgetl(f_id);
            obj.hor_L         = fscanf(f_id, '%f '); temp = fgetl(f_id);
            obj.hor_ShoS      = fscanf(f_id, '%f '); temp = fgetl(f_id);
            obj.hor_VolRat    = fscanf(f_id, '%f '); temp = fgetl(f_id);
            obj.hor_cel       = fscanf(f_id, '%f '); temp = fgetl(f_id);
            obj.hor_ceoc      = fscanf(f_id, '%f '); temp = fgetl(f_id);
            
            %- Dimensioni Piano Verticale
            temp = fgetl(f_id);
            tag_vert = 'Vertical:	';
            if ~strcmp(temp,tag_vert)
                error('Expected Vertical fields');
            end
            obj.vert.croot     = fscanf(f_id, '%f '); temp = fgetl(f_id); %croot
            obj.vert.croot_exp = fscanf(f_id, '%f '); temp = fgetl(f_id); %croot_exp
            obj.vert.ckink     = fscanf(f_id, '%f '); temp = fgetl(f_id); %ckink
            obj.vert.ctip      = fscanf(f_id, '%f '); temp = fgetl(f_id); %ctip
            obj.vert.x_apex    = nan;
            obj.vert.xle_exp   = fscanf(f_id, '%f '); temp = fgetl(f_id); %xle_exp
            obj.vert.xle_kink  = fscanf(f_id, '%f '); temp = fgetl(f_id); %xle_kink
            obj.vert.xle_tip   = fscanf(f_id, '%f '); temp = fgetl(f_id); %xle_tip
            obj.vert.y_root    = 0;
            obj.vert.y_rootexp = fscanf(f_id, '%f '); temp = fgetl(f_id); %y_root_exp
            obj.vert.y_kink    = fscanf(f_id, '%f '); temp = fgetl(f_id); %y_kink
            obj.vert_b         = fscanf(f_id, '%f '); temp = fgetl(f_id);
            obj.vert.y_tip     = obj.vert_b * 0.5;
            obj.vert_S         = fscanf(f_id, '%f '); temp = fgetl(f_id);
            obj.vert_sweepLE   = fscanf(f_id, '%f '); temp = fgetl(f_id);
            obj.vert_sweep025  = fscanf(f_id, '%f '); temp = fgetl(f_id);
            obj.vert_dihedral  = fscanf(f_id, '%f '); temp = fgetl(f_id);
            obj.vert_AR        = fscanf(f_id, '%f '); temp = fgetl(f_id);
            obj.vert_TR        = fscanf(f_id, '%f '); temp = fgetl(f_id);
            obj.vert.toc_root  = fscanf(f_id, '%f '); temp = fgetl(f_id);
            obj.vert.toc_kink  = fscanf(f_id, '%f '); temp = fgetl(f_id);
            obj.vert.toc_tip   = fscanf(f_id, '%f '); temp = fgetl(f_id);
            obj.vert.toc_avg   = fscanf(f_id, '%f '); temp = fgetl(f_id); %sarebbe t/c_avg
            obj.vert_mac       = fscanf(f_id, '%f '); temp = fgetl(f_id);
            obj.vert_L         = fscanf(f_id, '%f '); temp = fgetl(f_id);
            obj.vert_ShoS      = fscanf(f_id, '%f '); temp = fgetl(f_id);
            obj.vert_VolRat    = fscanf(f_id, '%f '); temp = fgetl(f_id);
            obj.vert_crud      = fscanf(f_id, '%f '); temp = fgetl(f_id);
            obj.vert_croc      = fscanf(f_id, '%f '); temp = fgetl(f_id);
            
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
            obj.TO_ISA_ST_SL = fscanf(f_id, '%f '); temp = fgetl(f_id);
            obj.TO_ISA_15_SL  = fscanf(f_id, '%f '); temp = fgetl(f_id);
            obj.TO_ISA_ST_5000 = fscanf(f_id, '%f '); temp = fgetl(f_id);
            obj.TO_ISA_15_5000    = fscanf(f_id, '%f '); temp = fgetl(f_id);
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
            obj.T0 = fscanf(f_id, '%f '); temp = fgetl(f_id);
            obj.Tmax_cont = fscanf(f_id, '%f '); temp = fgetl(f_id);
            catch ERR
                warning( ['Reading from file failed ',newline,ERR.message] )
            end
            fclose(f_id);
            
            %% Calcoli Derivati
            %_Pesi
            obj = obj.weightEst;
            %- Geometria delle Superfici Portanti
            try
                [obj.wing,obj.wing_S,obj.wing_AR,obj.wing_TR,obj.wing_sweepLE,obj.wing_025sweep,obj.wing_mac] = ...
                    obj.surfGeom( obj.wing,obj.wing_b,obj.wing_S,obj.wing_AR,obj.wing_TR,obj.wing_sweepLE );
            catch ERR_1
                warning( ' Error Calculating Wing Geomtery ' )
            end
            
            try
            [obj.hor,obj.hor_S,obj.hor_AR,obj.hor_TR,obj.hor_sweepLE,obj.hor_sweep025,obj.hor_mac] = ...
                obj.surfGeom( obj.hor,obj.hor_b,obj.hor_S,obj.hor_AR,obj.hor_TR,obj.hor_sweepLE );
            catch ERR_2
                warning( ' Error Calculating Horizontal Geomtery ' )
            end
            
            try
            [obj.vert,obj.vert_S,obj.vert_AR,obj.vert_TR,obj.vert_sweepLE,obj.vert_sweep025,obj.vert_mac] = ...
                obj.surfGeom( obj.vert,obj.vert_b,obj.vert_S,obj.vert_AR,obj.vert_TR,obj.vert_sweepLE );
            catch ERR_3
                warning( ' Error Calculating Verical Geomtery ' )
            end
            
            %Estrazione Aerodinamica
            obj.W_cr = 0.992*0.996*0.996*0.990*obj.MTOM*9.81; % Peso inizio crociera
            [~, ~, ~, obj.rho_cr] = atmosisa(obj.h_cr); % Densità inizio atmosfera
            %[obj.Cd0_cr,obj.e_cr] = obj.aero_finda;
            
            % Atro
            obj.Mark = obj.marker_assign();
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

                % Calcolo di Massa di Carburante e Zero Fuel Mass

                if isempty( obj.MPM )&& isempty( obj.OEM )
                    warining('Unable to calculate Operative Empty Mass');
                    obj.OEM = nan;
                else
                    if isempty( obj.MPM )
                        obj.MPM = obj.MZF - obj.OEM;
                    else
                        check_temp = obj.MZF - obj.MPM;
                        if isempty( obj.OEM )
                            obj.OEM = check_temp;
                        else
                            obj.checkLog( obj.OEM,check_temp,'OEM' );
                        end
                    end
                end
                
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
        
        function [lsurf_str,S,AR,TR,sweepLE,sweep025,cmac] = surfGeom(obj,lsurf_str,b,S,AR,TR,sweepLE_in)
            %surfGeom: funzione che calcola i parametri geometrici di una
            %superficie portante
            % INPUT:
            %   lsufr_str: struct contenente come campi i parametri
            %       geometrici della superficie%    
            %   b: aperura alare
            %   S: superficie alare
            TRi = 0;
            if isempty( lsurf_str.ckink )
                TRi = 1; lsurf_str.y_kink = 0;
            else
                if isempty( lsurf_str.y_kink )
                    warning( 'Kink Position Required, cannot proceed' );
                    lsurf_str.y_kink = nan;
                end
            end
            chk  = TRi == 1;
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
        
%         function [Cd0,e] = aero_finda(obj)
%             x0(1)    = 0.9; x0(2) = 0.00180;
%             x0(3:10) = [ obj.MTOM,obj.T0,obj.rho_cr, 0.75, obj.wing_b, obj.wing_AR,obj.wing_S,obj.V_cr ];
%             nROCs    = length(obj.ROC.roc(:));
%             x0(11)   = nROCs;
%             x0(12:11+nROCs) = obj.ROC.roc(:)';  x0(11+nROCs+1:11+2*nROCs) = obj.ROC.h(:)';
%             options = optimoptions('fsolve','algorithm','levenberg-marquardt');
%             x        = fsolve(@cruise_eq_sys_fun,x0,options);
%             Cd0      = x(2); e = x(1);
%         end
        
%         function y = cruise_eq_sys(x)
%             % x(1) = e
%             % x(2) = Cd0
%             % x(3) = MTOM,Vcr,T0,rho_cr,throttle_cr,b,AR,S,V_cr
%             MTOM = x(3);  T0 = x(4); rho_cr = x(5); throttle_cr = x(6);
%             b = x(7); AR = x(8); S = x(9); V_cr = x(10);
%             
%             sig = rho_cr/1.225; % rho_cr/rho_0
%             W_cr = 0.992*0.996*0.996*0.990*MTOM*9.81; % cruise weight [N]
%             T_cr = T0*sig*throttle_cr; % Spinta in crociera [N]
%             CL   = obj.W_cr*2/( rho_cr*V_cr^2*S ); % CL di cruise
%             K    = 1/( pi*AR*x(1) ); f = x(2)*S; be = b*sqrt(x(1));
%             y(1) = 0.5*rho_cr*V_cr*S*( x(2) + K*CL^2 ) - T_cr; % T = D
%             y(2) = 1.54*(T_cr/W_cr)*sqrt(T_cr/f)/sig - 2.2*(W_cr/9.81)/(be^2) / ( sqrt( T_cr/9.81 * sig/f ) ); % RC_max -> CONTROLLA
%         end
        
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
