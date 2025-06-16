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
        yob_kink     % Positions along span
        eps          % Twist Distrubution [deg]
        toc          % t/c distribution
        c            % chord distribution [m]
        %% Initial Guess Aerodynamic Data
        CL_cr        % CL cruise estimated
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
            c = 0.0199; d = 0.7531;

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
            /( V_cr_vet( ch(1) )^2*rho );
        end
        
        function fig = preliminary_polar_plot( obj )
           fig    = figure('Name','Preliminary Polar Estimation');
           ax_1   = subplot(3,2,1,'Parent',fig); hold( ax_1,'on' );
           colors = [
               1, 0, 0;    % Rosso
               0, 1, 0;    % Verde
               0, 0, 1;    % Blu
               1, 1, 0;    % Giallo
               0, 1, 1;    % Ciano
               1, 0, 1;    % Magenta
               1, 0.647, 0 % Arancione
               ];
           [T, a_sound, P, rho] = atmosisa(obj.TLARs.cruise.h); M_crit = obj.TLARs.cruise.M;
           %% Aircraft Polar
           CL_vet   = linspace( -0.1*obj.CLmax_cr,obj.CLmax_cr,500); CL_vet = CL_vet(:);
           CD_vet   = obj.SizHis(end-1).CD0 + CL_vet.^2 /( pi*obj.ARw*obj.TLARs.e);
           lin(1,1) = plot( ax_1,CD_vet,CL_vet ); lin(1,1).LineStyle = '-'; lin(1,1).LineWidth = 2;
           lin(1,1).Annotation.LegendInformation.IconDisplayStyle = 'off';
           npts = 5; PT = {'E','P','A','stall','cruise';'o','o','o','o','square'}; PTO(1).CL   = sqrt( pi*obj.ARw*obj.TLARs.e*obj.SizHis(end-1).CD0 ); PTO(1).CD = 2*obj.SizHis(end-1).CD0; % Pto E
           i = 2; PTO(i).CL = sqrt( 3*pi*obj.ARw*obj.TLARs.e*obj.SizHis(end-1).CD0 ); PTO(i).CD = 4*obj.SizHis(end-1).CD0; % Pto P
           i = 3; PTO(i).CL = sqrt( 1/3*pi*obj.ARw*obj.TLARs.e*obj.SizHis(end-1).CD0 ); PTO(i).CD = 4/3*obj.SizHis(end-1).CD0; % Pto P
           i = 4; PTO(i).CL = obj.CLmax_cr; PTO(i).CD = obj.SizHis(end-1).CD0 + obj.CLmax_cr.^2 /( pi*obj.ARw*obj.TLARs.e);
           i = 5; PTO(i).CL = obj.CL_cr; PTO(i).CD = obj.SizHis(end-1).CD0 + obj.CL_cr.^2 /( pi*obj.ARw*obj.TLARs.e);
           for i = 1:npts
               lin(1,i+1) = plot( ax_1,PTO(i).CD,PTO(i).CL ); lin(1,i+1).LineStyle = 'none'; lin(1,i+1).Marker = PT{2,i};
               lin(1,i+1).MarkerSize = 6; lin(1,i+1).MarkerEdgeColor = colors(i,:); lin(1,i+1).LineWidth = 1.2;
               lin(1,i+1).DisplayName = ['Point ',PT{1,i}];
           end
           legend( ax_1,'Interpreter','Latex'); title('Preliminary Aircraft Polar in Cruise','Interpreter','Latex'); xlabel( 'C$_D$','Interpreter','Latex' ); ylabel( 'C$_L$','Interpreter','Latex' ); 
           %% CL - V Plot
           % Cruise at max height and intermediate weight
           ax_2     = subplot(3,2,3:4,'Parent',fig); j = 2; i = 1;
           V_cr     =  sqrt( 9.81*obj.SizHis(end-1).WoS*0.5*( obj.MCroMTo(1)+obj.MCroMTo(2) )*2 ./ (rho.*CL_vet( CL_vet>0.1 ) ) ) ;
           lin(2,1) = plot( ax_2,V_cr,CL_vet( CL_vet>0.1 ) ); hold( ax_2,'on' );
           lin(j,i).LineStyle = '-'; lin(j,i).LineWidth = 2; lin(j,i).Annotation.LegendInformation.IconDisplayStyle = 'off';
           ivt = 1:npts;
           for i = 1:4
               PTO(i).V  = sqrt( 2*9.81*obj.SizHis(end-1).WoS*0.5*( obj.MCroMTo(1)+obj.MCroMTo(2) )/( PTO(i).CL*rho ) );
               lin(j,i)  = plot( ax_2,PTO(i).V,PTO(i).CL ); lin(j,i).LineStyle = 'none'; lin(j,i).Marker = PT{2,i};
               lin(j,i).MarkerSize  = 6; lin(j,i).MarkerEdgeColor = colors(i,:); lin(j,i).LineWidth = 1.2;
               lin(j,i).DisplayName = ['Point ',PT{1,i}];
           end
           i = 5; PTO(i).V  = obj.TLARs.cruise.M*a_sound; lin(j,i)  = plot( ax_2,PTO(i).V,PTO(i).CL ); lin(j,i).LineStyle = 'none'; lin(j,i).Marker = PT{2,i};
           lin(j,i).MarkerSize  = 6; lin(j,i).MarkerEdgeColor = colors(i,:); lin(j,i).LineWidth = 1.2;
           lin(j,i).DisplayName = ['Point ',PT{1,i}]; xlabel( 'V [m/s] ','Interpreter','Latex' ); ylabel( 'C$_L$','Interpreter','Latex' ); 
           legend( ax_2,'Interpreter','Latex'); title('V - C$_L$ Diagram','Interpreter','Latex');
           lin(j,i).LineStyle = 'none'; lin(j,i).Marker = 'square';
           %% Preq - V
           % Cruise at max height and intermediate weight
           ax_3 = subplot(3,2,5:6,'Parent',fig); j = 3; i = 1;
           D    = 0.5*rho*V_cr.^2*obj.Sw.*( obj.SizHis(end-1).CD0 + ...
               1/(pi*obj.ARw*obj.TLARs.e)*( 9.81*obj.SizHis(end-1).WoS*0.5*( obj.MCroMTo(1)+obj.MCroMTo(2) )*2 ./ (rho.*V_cr.^2 ) ).^2  + ...
               20*( V_cr./a_sound - M_crit ) ).*( V_cr./a_sound > M_crit );
           Preq = D.*V_cr; lin(j,i) = plot( ax_3,V_cr,Preq ); hold( ax_3,'on' );
           lin(j,i).LineStyle = '-'; lin(j,i).LineWidth = 2; lin(j,i).Annotation.LegendInformation.IconDisplayStyle = 'off';
           for i = 1:5
               PTO(i).D  = 0.5*rho*PTO(i).V.^2*obj.Sw.*( obj.SizHis(end-1).CD0 + ...
               1/(pi*obj.ARw*obj.TLARs.e)*( 9.81*obj.SizHis(end-1).WoS*0.5*( obj.MCroMTo(1)+obj.MCroMTo(2) )*2 ./ (rho.*PTO(i).V.^2 ) ).^2 );
               lin(j,i)  = plot( ax_3,PTO(i).V,PTO(i).D*PTO(i).V ); lin(j,i).LineStyle = 'none'; lin(j,i).Marker = PT{2,i};
               lin(j,i).MarkerSize  = 6; lin(j,i).MarkerEdgeColor = colors(i,:); lin(j,i).LineWidth = 1.2;
               lin(j,i).DisplayName = ['Point ',PT{1,i}];
           end
           legend( ax_3,'Interpreter','Latex'); title('V - $\Pi_{req}$ Diagram','Interpreter','Latex');  xlabel( 'V [m/s] ','Interpreter','Latex' ); ylabel( '$\Pi_{req}$ [W]','Interpreter','Latex' ); 
           lin(j,i).LineStyle = 'none'; lin(j,i).Marker = 'square';
           %% E - V
           ax_4 = subplot(3,2,2,'Parent',fig); j = 4; i = 1;
           E_v =  9.81*obj.SizHis(end-1).WoS*0.5*( obj.MCroMTo(1)+obj.MCroMTo(2) )*2 ./ ( rho*V_cr.^2 )... % CL
            ./ ( obj.SizHis(end-1).CD0 + ...
               1/(pi*obj.ARw*obj.TLARs.e)*( 9.81*obj.SizHis(end-1).WoS*0.5*( obj.MCroMTo(1)+obj.MCroMTo(2) )*2 ./ (rho.*V_cr.^2 ) ).^2 );
            lin(j,i) = plot( ax_4,V_cr,E_v ); hold( ax_4,'on' );
           lin(j,i).LineStyle = '-'; lin(j,i).LineWidth = 2; lin(j,i).Annotation.LegendInformation.IconDisplayStyle = 'off';
           for i = 1:5
               lin(j,i)  = plot( ax_4,PTO(i).V,PTO(i).CL/PTO(i).CD ); lin(j,i).LineStyle = 'none'; lin(j,i).Marker = PT{2,i};
               lin(j,i).MarkerSize  = 6; lin(j,i).MarkerEdgeColor = colors(i,:); lin(j,i).LineWidth = 1.2;
               lin(j,i).DisplayName = ['Point ',PT{1,i}];
           end
           legend( ax_4,'Interpreter','Latex'); title('V - E Diagram','Interpreter','Latex');  xlabel( 'V [m/s] ','Interpreter','Latex' ); ylabel( 'E','Interpreter','Latex' ); 
           lin(j,i).LineStyle = 'none'; lin(j,i).Marker = 'square';
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
        function CDlow = CDlow_Mach(obj,alpha,alpha_v,cds)
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

            cdfit = polyval( obj.low_speed.meanprofile.poly_drag,alpha ); % cd del profilo 2d ottenuto con il polinomio
            CL    = obj.low_speed.lift_eval( alpha,obj.low_speed.prf3DClean );
            u     = obj.interpolateFromCSV('AR', obj.TRw, obj.ARw);
            v     = obj.interpolateFromCSV('TR', obj.ARw, obj.TRw);
            w     = obj.interpolateFromCSV('ctcr', obj.ARw, obj.TRw);


            t1 = CL.^2./(pi*obj.ARw*u);%c'è un fattore s che non sappiamo cosa significa anche nell'excel non viene calcolato
            t2 = v*CL*obj.low_speed.meanprofile.eps_ae*obj.low_speed.meanprofile.a;
            t3 = ( obj.low_speed.meanprofile.eps_ae*obj.low_speed.meanprofile.a )^2*w; % eps_ae + Cla

            CDi = t1(:)+t2(:)+t3(:);

            CDlow = cdfit(:) + CDi(:); %non c'è contributo di wave perchè siamo a basso mach

        end
        
        function value = interpolateFromCSV(obj,name,xq, yq)
            %valuta u, v, w
            %INput:
            %pattern = l'inizio del nome del file csv: può essere AR*(u),TR*(v),ctcr*(w)
            %xq,yq = punto query in cui voglio u,v,w
            files = [ obj.main_fold,'\functions\Wing_Design\u_TR\',name,'_file.csv'];
            MAT = readmatrix(files);
            F = scatteredInterpolant( MAT(:,1),MAT(:,3),MAT(:,2) );
            value = F(xq,yq);
        end
        
        function CDfun = CDtransonic( obj, alpha, M_cruise, alpha_v, cds )
            % CDw Calcola il coefficiente di drag totale di un'ala in cruise
            %
            % INPUT:
            %   alpha       - vector of alpha
            %   M_cruise    - cruise mach
            %
            % OUTPUT:
            %   CDfun       - Drag totale
            ka = 0.95; % TEMPORARYYYYY
            sweep_c4 = obj.low_speed.sweep*pi/180; tc_mean = obj.low_speed.meanprofile.tc;
            if obj.high_speed.meanprofile.M ~= M_cruise
                warning('M used for calculations is different from Mach at which the aerodynamic data are calculated');
            end
            MDD = ka/(cos(sweep_c4)) - ( tc_mean / ( cos(sweep_c4)^2) )-...
                ( obj.CL_cr/(10*(cos(sweep_c4))^3 ) );
            Mcrit = MDD-(0.1/80)^(1/3);

            DeltaCd_wave = 20*( M_cruise-Mcrit )^4;

            %calcolo dei coefficienti u,v,z per il CDi con la function
            %interpolatefromcsv fatta a parte
            u = obj.interpolateFromCSV('AR',   obj.TRw, obj.ARw);
            v = obj.interpolateFromCSV('TR',   obj.ARw, obj.TRw);
            w = obj.interpolateFromCSV('ctcr', obj.ARw, obj.TRw);
            
            % Wing Cd, CL
            
            
            if isnan( obj.high_speed.meanprofile.poly_drag )
                % If poly_drag object is not initialized, it initializes
                % it.
                obj.high_speed.meanprofile.poly_drag = obj.high_speed.poly_drag( alpha_v,cds );
            end
            CLw_c = obj.high_speed.lift_eval( alpha,obj.high_speed.prf3DClean ); 
            Cd_avg = polyval( obj.high_speed.meanprofile.poly_drag,alpha );

            t1 = CLw_c.^2 ./ (pi*obj.ARw*u); % c'è un fattore s che non sappiamo cosa significa anche nell'excel non viene calcolato
            t2 = v*CLw_c*obj.high_speed.meanprofile.eps_ae*obj.high_speed.meanprofile.a; % Cl_alpha_m;
            t3 = ( obj.high_speed.meanprofile.eps_ae*obj.high_speed.meanprofile.a )^2*w;

            CDi   = t1+t2+t3;

            CDfun = Cd_avg(:) + CDi(:) + DeltaCd_wave(:);
        end

        function Mdd(obj)
            % Mdd_fun Fa il check sulla buffet barrier
            %
            % INPUT:
            %
            % OUTPUT:
            %
            v1 = linspace(0,0.5,11);
            v2 = linspace(0.52,0.87,9);
            M_vett = [v1,v2]; % Vettore di mach assunto
            Cl_max = obj.low_speed.prf3DClean.clmax *0.8;% prende il Clmax 2d e lo moltiplica per 0.8 per ottenere quello 3d
            tc_mean = obj.low_speed.meanprofile.tc;
            cosc4 = cos(obj.low_speed.sweep/57.3); %cos dell'angolo di freccia a c/4
            deltamcc = 0.06; % delta mach critico per profili supercritici(0.06)
            Cl_d0 = Cl_max./sqrt(1-M_vett.^2);% CL con correzione di prandtl-glauert

            x = (tc_mean/cosc4);
            K1 = 2.8355*x^2-1.9072*x+0.9499;
            K2 = 0.2*(1-2.131*x);
            Cl_mcc = (K1*cosc4^2-(M_vett-deltamcc)*cosc4^3)/K2;
            diff = Cl_mcc -Cl_d0;
            jstart = find(diff<0.01,1,"first");
            Msub = M_vett(jstart:end);
            for j = 1:length(Msub)
                Mdd(j) = (Msub(j)-0.06)/(1.02+0.08*(1-cosc4));
                Cl_mdd(j) = (K1*cosc4^2-Mdd(j)*cosc4^3)/K2;
            end
            fig = figure( 'Name','Buffet Check' ); ax_b = axes('Parent',fig);
            plot( ax_b,M_vett(1:jstart),Cl_d0(1:jstart)); hold(ax_b,'on');
            plot( ax_b,Msub,Cl_mcc(jstart:end) );
            plot( ax_b,Msub,Cl_mdd );
            plot( ax_b,obj.TLARs.cruise.M,obj.CL_cr,'o' ); %current point
            legend('Cl a M diverso da 0','Cl a Mcc','Cl a MDD','Current point')
            % Chiedi all'utente se vuole continuare
            risposta = input('Do you want to go on? (s/n): ', 's');
            % Controlla la risposta
            if risposta == 's'
                disp('It will go on. ');
            elseif risposta == 'n'
                disp('It will stop.');
            else
                disp('No valid answer.');
            end
        end

        function plot_fun(obj,alfavett,CL_low,CL_high,CD_low,CD_high)
            %plot effettua i grafici mettendo a confronto CL-alfa, CD-CL
            %
            %INPUT:
            %
            % alphavett = vettore di alfa
            % CL_low = CL a basse velocità

            
            fig = figure("Name","Wing Polars");
            X = {alfavett,alfavett,CL_low,CL_high};
            Y = {CL_low,CL_high,CD_low,CD_high};
            
            for i = 1:4 %1 a 6 se metti CM
                ax_c(i) = subplot(2, 2, i,"Parent",fig); 
                plot(X{i}, Y{i}, 'LineWidth', 1.5);
            end

        end
    end
end

