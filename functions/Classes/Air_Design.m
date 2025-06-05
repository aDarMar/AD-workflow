classdef Air_Design
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        main_fold
        % Struct with TLARs
        TLARs
        %% Weights and Weight Franctions
        % Class 0 Weight Estimation
        MTOM_est0   %           [Kg]
        EM_est0     % Empy Mass [Kg]
        Mpay_max    % Max Payload Mass [Kg]
        M_crew      % Crew Mass [Kg]
        % Fuel Fractions
        Mff         % W_fin/W_in for all the mission
        MLndoMTo    % Landing Mass / Take-Off Mass
        MCroMTo     % Take-Off Mass / Cruise Mass [initial,final,avg]
        % Linear MTOMvsME cruve coefficients
        weight_coeff_c
        weight_coeff_d
        %% Initial Guess Geometric Parameters
        SizHis      % Sizing History Vector 
        dCD0_wave
        T0          % Initial Guess Max Static Thrust [N] 
        % First-Guess Geometric Parameters
        Sw
        bw
        ARw          % Assumed Aspect-Ratio
        TRiw         % Assumed Taper C_kink/C_root
        TRw          % Taper C_tip/C_root
        sweepw       % Assumed L.E. Sweep [deg]
        dihedralw    % Assumed Dihedral Angle [deg]
        AioSw        % Ratio Between inner panel Area over Wing Area
        iw           % Wing Incidence Angle [deg]
        wingapex     % Wing Apex Coordinates [x,y,z] [m]
        % [root,kink,tip]
        yob_kink      % Positions along span
        eps         % Twist Distrubution [deg]
        toc         % t/c distribution
        c           % chord distribution [m]
        %% Initial Guess Aerodynamic Data
        CL_cr       % CL cruise estimated
        CLmax_cr
        CLmax_TO
        CLmax_LND

        equiv_wing  % Wing class object representing the equivalent wing
        low_speed
        high_speed
    end
    
    methods
        function obj = Air_Design(name_file,main_path)
            %UNTITLED Construct an instance of this class
            %   name_file: path of TLARs text file
            obj.main_fold = main_path;
            obj.dCD0_wave    = 0.0015;
            %% Reading TLARS
            obj = obj.read_TLARs( name_file );
            obj = obj.read_geometry;
            obj.TLARs.AR = obj.ARw;
            %% Fuel Fractions
            % Add choice for Mff type (Breguet and statistical) for climb
            [obj.Mff,temp2,obj.MLndoMTo,obj.MCroMTo] = obj.fuel_fraction;
            %% Weight Estimation
            obj.Mpay_max = obj.TLARs.npax*215/2.2046;                       %[Kg]
            obj.M_crew   = ( obj.TLARs.ncrew+obj.TLARs.npil )*205/2.2046;   %[Kg]
            Mres = 0; Mfo= 0;
            obj.weight_coeff_c = 1 - (1+Mres)*(1-obj.Mff) - Mfo; 
            obj.weight_coeff_d = obj.Mpay_max + obj.M_crew;
        end
        %% Reading Input Files
        function obj = read_TLARs(obj,tlars_file_path)
            %read_TLARs function that reads TLARs and assumptions for a design airplane
            %from a TXT file
            %   tlars_file_path: string containing TLARs txt file path
            %   TLARS: struct containing the TLARs
            %
            %   02/05/25 v. 0.1e

            f_id = fopen(tlars_file_path,'r');
            % Grafica di Lettura
            disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
            disp('Reading TLARs ...' )
            disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
            try
                temp = fgetl(f_id); tag = 'TLARS';
                if ~strcmp(temp,tag)
                    error(' The file opened does not contain TLARs')
                end
                % General TLARs
                temp = fgetl(f_id); tag = 'General';
                if ~strcmp(temp,tag)
                    error(' Expected General TLARS')
                end

                obj.TLARs.npax  = fscanf(f_id,'%f '); temp = fgetl(f_id);
                obj.TLARs.ncrew   = fscanf(f_id,'%f ');
                if isempty( obj.TLARs.ncrew ) || obj.TLARs.ncrew < ceil( obj.TLARs.npax/50 )
                    obj.TLARs.ncrew = ceil( obj.TLARs.npax/50 );
                end
                temp = fgetl(f_id);
                obj.TLARs.npil    = fscanf(f_id,'%f '); temp = fgetl(f_id);

                temp = fgetl(f_id); tag = 'Aerodynamics';
                if ~strcmp(temp,tag)
                    error(' Expected Aerodynamics TLARs')
                end
                obj.TLARs.e          = fscanf(f_id,'%f '); temp = fgetl(f_id);
                obj.TLARs.dCD0_f_TO  = fscanf(f_id,'%f '); temp = fgetl(f_id);
                obj.TLARs.dCD0_f_LND = fscanf(f_id,'%f '); temp = fgetl(f_id);
                temp  = fscanf(f_id,'%f ');
                if isempty( temp )
                    temp = ( obj.TLARs.dCD0_f_TO+obj.TLARs.dCD0_f_LND )*0.5;
                end
                obj.TLARs.dCD0_f_App = temp; temp = fgetl(f_id);
                obj.TLARs.dCD0_lgs   = fscanf(f_id,'%f ');  temp = fgetl(f_id);
                obj.TLARs.de_TO      = -fscanf(f_id,'%f '); temp = fgetl(f_id);
                obj.TLARs.de_LND     = -fscanf(f_id,'%f '); temp = fgetl(f_id);

                temp = fgetl(f_id); tag = 'Propulsion';
                if ~strcmp(temp,tag)
                    error(' Expected Propulsive TLARS')
                end
                obj.TLARs.T0oTmc  = fscanf(f_id,'%f '); temp = fgetl(f_id);
                obj.TLARs.nengine = fscanf(f_id,'%f '); temp = fgetl(f_id);

                %% Take-Off
                temp = fgetl(f_id); tag = 'Take-Off';
                if ~strcmp(temp,tag)
                    error(' Expected Take-Off Requirements')
                end
                temp = fscanf(f_id,'%f '); obj.TLARs.TO.fieldmax = convlength(temp,'ft','m');
                temp = fgetl(f_id);

                %% Climb Specs
                temp = fgetl(f_id); tag = 'Climb';
                if ~strcmp(temp,tag)
                    error(' Expected Climb Requirements')
                end
                temp = fscanf(f_id,'%f '); obj.TLARs.climb.RoC = convvel(temp,'ft/min','m/s');
                temp = fgetl(f_id);
                temp = fscanf(f_id,'%f '); obj.TLARs.climb.V = convvel(temp,'kts','m/s');
                temp = fgetl(f_id);
                obj.TLARs.climb.E   = fscanf(f_id,'%f '); temp = fgetl(f_id);
                temp  = fscanf(f_id,'%f '); obj.TLARs.climb.cj = temp/3600;
                temp= fgetl(f_id);

                %% Cruise Specs.
                temp = fgetl(f_id); tag = 'Cruise';
                if ~strcmp(temp,tag)
                    error(' Expected Cruise Requirements')
                end
                temp = fscanf(f_id,'%f '); obj.TLARs.cruise.h = convlength(temp,'ft','m');
                temp = fgetl(f_id);
                temp = fscanf(f_id,'%f '); obj.TLARs.cruise.R = convlength(temp,'naut mi','m');
                temp = fgetl(f_id);
                [T, a_sound, P, rho] = atmosisa(obj.TLARs.cruise.h);
                obj.TLARs.cruise.M = fscanf(f_id,'%f '); temp = fgetl(f_id);
                obj.TLARs.cruise.V = obj.TLARs.cruise.M*a_sound;
                obj.TLARs.cruise.E = fscanf(f_id,'%f '); temp = fgetl(f_id);
                temp  = fscanf(f_id,'%f '); obj.TLARs.cruise.cj = temp/3600;
                temp = fgetl(f_id);

                %% Loiter Specs.
                temp = fgetl(f_id); tag = 'Loiter';
                if ~strcmp(temp,tag)
                    warning('Loiter not Defined')
                    obj.TLARs.loiter.End = 0;
                else
                    obj.TLARs.loiter.E = fscanf(f_id,'%f '); temp = fgetl(f_id);
                    temp = fscanf(f_id,'%f '); obj.TLARs.loiter.cj = temp/3600;
                    temp = fgetl(f_id);
                    temp = fscanf(f_id,'%f '); obj.TLARs.loiter.End = temp*60;
                    temp = fgetl(f_id);
                end
                %% Alternate Specs
                temp = fgetl(f_id); tag = 'Alternate';
                if ~strcmp(temp,tag)
                    warning('Alternate not Defined')
                    obj.TLARs.alter.R = 0;

                    obj.TLARs.alter.h = 1; obj.TLARs.alter.M = 1; % continua


                else
                    temp = fscanf(f_id,'%f '); obj.TLARs.alter.h = convlength(temp,'ft','m');
                    temp = fgetl(f_id);
                    temp = fscanf(f_id,'%f '); obj.TLARs.alter.R = convlength(temp,'naut mi','m');
                    temp = fgetl(f_id);
                    [T, a_sound, P, rho] = atmosisa(obj.TLARs.alter.h);
                    obj.TLARs.alter.M = fscanf(f_id,'%f '); temp = fgetl(f_id);
                    obj.TLARs.alter.V = obj.TLARs.alter.M*a_sound;
                    obj.TLARs.alter.E = fscanf(f_id,'%f '); temp = fgetl(f_id);
                    temp  = fscanf(f_id,'%f '); obj.TLARs.alter.cj = temp/3600;
                    temp = fgetl(f_id);

                end

                %% Landing Specs
                temp = fgetl(f_id); tag = 'Landing';
                if ~strcmp(temp,tag)
                    error('Expected Landing TLARS')
                end
                temp = fscanf(f_id,'%f '); obj.TLARs.LND.SGmax = convlength(temp,'ft','m');
                temp = fgetl(f_id);

                fclose(f_id);
                disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
            catch ERR
                fclose('all');
                error( ['An error occurred while reading obj.TLARs, cannot proceed ...',ERR.message] );
            end


        end
        
        function obj = read_geometry(obj)
            pth = [obj.main_fold,'\tlars\Geometry.txt'];
            f_id = fopen(pth);
            % First Line Check
            temp = fgetl(f_id);
            tag = 'Design Parameters';
            if ~strcmp(temp,tag)
                error('File format inavalid')
            end
            
            % Wing Check
            temp = fgetl(f_id);
            tag = 'Wing';
            if ~strcmp(temp,tag)
                error('Expected Wing Fields')
            end
            temp = fgetl(f_id);
            tag = 'General';
            if ~strcmp(temp,tag)
                error('Expected Wing Fields')
            end

            obj.ARw    = fscanf(f_id,'%f '); temp = fgetl(f_id);
            obj.TRw    = fscanf(f_id,'%f '); temp = fgetl(f_id);
            obj.TRiw   = fscanf(f_id,'%f '); temp = fgetl(f_id);
            obj.AioSw   = fscanf(f_id,'%f '); temp = fgetl(f_id);
            obj.sweepw    = fscanf(f_id,'%f '); temp = fgetl(f_id);
            obj.dihedralw = fscanf(f_id,'%f '); fgetl(f_id);
            obj.iw = fscanf(f_id,'%f '); fgetl(f_id);
            temp = fscanf(f_id,'%f '); fgetl(f_id);
            obj.wingapex.x = temp(1);
            obj.wingapex.y = temp(2);
            obj.wingapex.z = temp(3);
            % temp = fgetl(f_id);
            % tag = 'Sections';
            % if ~strcmp(temp,tag)
            %     error('Expected Wing Fields')
            % end
            % temp = fgetl(f_id);
            % 
            % obj.eps = fscanf(f_id,'%f '); % gets the first row
            % %temp = textscan(f_id, '%f', 'Delimiter', '\t');
            fclose(f_id);
        end
        %% Class 0 Weight Estimation
        function [Mff_stat,Mff_breg,MLNDoMTO,WTOoWcr] = fuel_fraction(obj)
            %UNTITLED2 Summary of this funcjion goes here
            %   1 - engine start and warm-up
            %   2 - taxi
            %   3 - take-off
            %   4 - climb (statistical)
            %   5 - cruise (for given range)
            %   6 - loiter (in terms of endurance)
            %   7 - descend
            %   8 - Landing, Taxi and Shutdown
            %   9 - Alternate (for a given range)
            % Definizione del Profilo di Volo
            Wfrac = nan(9,1);

            Wfrac(1) = 0.99;
            Wfrac(2) = 0.99;
            Wfrac(3)= 0.995;
            Wfrac(4) = 0.980;
            Wfrac(7) = 0.99;
            Wfrac(8) = 0.992;

            % Breguet
            Wfrac(5) = exp( - obj.TLARs.cruise.R*obj.TLARs.cruise.cj/( obj.TLARs.cruise.V * obj.TLARs.cruise.E ) ); % Cruise
            Wfrac(6) = exp( - obj.TLARs.loiter.cj*obj.TLARs.loiter.End/obj.TLARs.loiter.E ); % Loiter
            Wfrac(9) = exp( - obj.TLARs.alter.R*obj.TLARs.alter.cj/( obj.TLARs.alter.V * obj.TLARs.alter.E ) ); % Alternate
            n_phases = 9;

            Mff_stat = 1;
            for i = 1:n_phases
                Mff_stat = Mff_stat*Wfrac(i);
            end

            obj.TLARs.climb.End = obj.TLARs.cruise.h/obj.TLARs.climb.RoC; obj.TLARs.climb.R = obj.TLARs.climb.End*sqrt( obj.TLARs.climb.V^2 - obj.TLARs.climb.RoC^2 );
            Wfrac(4) = exp( - obj.TLARs.climb.cj*obj.TLARs.climb.End/obj.TLARs.climb.E ); % Climb

            Mff_breg = 1;
            for i = 1:n_phases
                Mff_breg = Mff_breg*Wfrac(i);
            end
            MLNDoMTO = Mff_breg/Wfrac(8);
            WTOoWcr(1) = Wfrac(1)*Wfrac(2)*Wfrac(3)*Wfrac(4);
            WTOoWcr(2) = Wfrac(1)*Wfrac(2)*Wfrac(3)*Wfrac(4)*Wfrac(5);
        end
        %% Sizing Point
        function [Cd0,Swet] = polar_est(obj,S,MTOM)
            %polar_est Function that evaluates the polar with a statistical approach.
            %   MTOM:   in Kg
            %   S   :   estimated wing surface [m]

            lb2kg = 0.45359237; ft2m = 0.3048;
            % From Roskam, values for Transport Jets
            % log10( Swet ) = c + d*log10( WTO ) with WTO in [lb], Swet in [ft^2]
            c = 0.0199; d = 0.7331;

            Swet  = 10^( c+d*log10(MTOM/lb2kg) ); %Swet in ft^2
            temp = readmatrix([obj.main_fold,'\statistical_data\cf_vs_Swet.csv']);
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

        function obj = final_out( obj,idxs,CLmax_TO_vett,CLmax_CR_vett,...
                CLmax_LND_vettiS,V_cr_vet,h_cr,iS )
            %final_out: saves the values chosen from the Sizing plot inside
            %the design wing object
            if nargin <8
                iS = length ( obj.SizHis(:) );
            end
            obj.Sw = obj.SizHis(iS).S;
            obj.bw = sqrt( obj.Sw*obj.ARw );
            obj.T0 = obj.SizHis(iS-1).ToW*obj.MTOM_est0*9.81;
            %% Aerodynamic Data
            i = 1; ch = idxs( idxs(:,1) == i,2 ); % selects the first choice made for CLmax@T/O
            obj.CLmax_TO  = CLmax_TO_vett( ch(1) );
            i = 2; ch = idxs( idxs(:,1) == i,2 ); % selects the first choice made for CLmax@LAND
            obj.CLmax_LND = CLmax_LND_vettiS( ch(1) );
            i = 3; ch = idxs( idxs(:,1) == i,2 ); % selects the first choice made for CLmax@cruise ( CL cruise is in the middle )
            obj.CLmax_cr  = CLmax_CR_vett( ch(2) );
            i = 4; ch = idxs( idxs(:,1) == i,2 );
            [T,a,P,rho]   = atmosisa( h_cr( ch(1) ) );
            obj.CL_cr     = 9.81*obj.SizHis(iS).WoS*...
                0.5*( obj.MCroMTo(1)+obj.MCroMTo(2) )*2 ...
            /( V_cr_vet(ch)^2*rho );
        end
        
        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
        %% Wing Design
        function obj = equivalent_wing_def(obj,ctip)
            %equivalent_wing_def: function that defines the equivalent wing
            %given as an object of WingClass class.
            yroot = 0;
            croot_eq = ((((obj.Sw/(0.5*obj.bw*(1-yroot/obj.bw)))-2*ctip)/(0.5*obj.bw*(1-yroot/obj.bw)))*0.5* ...
                obj.bw*yroot/obj.bw)+(obj.Sw/(0.5*obj.bw*(1-yroot/obj.bw)))-ctip;
            
            %sweep_eq    = atan((xLE_tip-xLE_root)/(ytip-yroot))*57.3;
            apexC = [0,0,0];
            
            sectsGeom = nan(2,8);
            i = 1; sectsGeom(:,i) = [croot_eq;ctip];
            sectsAero = nan(2,8);

            % slop_c      = (ctip_eq-croot_eq)*2/obj.bw;
            % slop_tw     = (eps_tip-eps_root)/ytip;
            % slop_xle    = (xLE_tip-xLE_root)/ytip;

            obj.equiv_wing = WingClass(obj.bw*0.5,obj.sweepw,obj.dihedralw,nan,apexC,nan,...
                sectsGeom, sectsAero); %Costruttore
        end
        % Preliminary Drag Estimation
        function CDlow = CDlow_Mach(obj,alpha,CL,alpha_v,cds)
            % CDSTALL Calcola il drag totale a basso Mach (senza drag d’onda)
            %
            % INPUT:
            %   alpha_v  - Vettore angoli d'attacco per il profilo 2D
            %   cds      - Vettore corrispondente dei cd 2D ordinati per
            %               [ cd(root,alpha_i),cd(kink,alpha_i),cd_kink,alpha_i) ]
            %   alpha    - Input alpha
            %   CL       - Portanza dell’ala in low Mach
            %
            % OUTPUT:
            %   CDstall  - Coefficiente di drag totale a low Mach
            
            if isnan( obj.low_speed.meanprofile.poly_drag )
                % If poly_drag object is not initialized, it initializes
                % it.
                obj.low_speed.meanprofile.poly_drag = obj.low_speed.poly_drag( alpha_v,cds );
            end

            % erroreMAX = 0.05; % imposto il minimo errore
            % for n=1:5
            %     p     = polyfit(a,cd_mean,n);
            cdfit = polyval( obj.low_speed.meanprofile.poly_drag,alpha ); % cd del profilo 2d ottenuto con il polinomio
            %     err   = norm(cd_mean-cdfit);
            %     if err<erroreMAX
            %         break;
            %     end
            % end
            u = interpolateFromCSV('AR*.csv', obj.TRw, obj.ARw);
            v = interpolateFromCSV('TR*.csv', obj.ARw, obj.TRw);
            w = interpolateFromCSV('ctcr*.csv', obj.ARw, obj.TRw);


            t1 = CL^2/(pi*obj.ARw*u);%c'è un fattore s che non sappiamo cosa significa anche nell'excel non viene calcolato
            t2 = v*CLw_s*obj.low_speed.meanprofile.eps_ae*obj.low_speed.meanprofile.a;
            t3 = ( obj.low_speed.meanprofile.eps_ae+obj.low_speed.meanprofile.a )^2*w; % eps_ae + Cla

            CDi = t1+t2+t3;

            CDlow = cdfit+ CDi; %non c'è contributo di wave perchè siamo a basso mach

        end

    end
end

