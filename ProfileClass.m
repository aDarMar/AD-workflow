classdef ProfileClass %<handle
    properties
        flag
        % Geometry
        c          %corda del profilo [m]
        eps                        %svergolamento [deg]
        tc           %Spessore medio
        LERc         % Raggio del L.E. adim
        xtc             % Coordinata adim di max spessore
        xtrUp        % Lunghezza del BL laminare su dorso/ventre
        xtrLow 
        dY                           % Altezza tra i punti a 0.005c e 0.1c
        %Aerodynamics
            % Lift
        h                            % Quote di Volo [m]
        M                           % Mach di volo a cui sono riferite le grandezze aerodinamiche
        a 
        cl0 
        clstar 
        clmax 
        alphamax 
        alpha0l 
        alphastar 
        cmac 
            % Drag: non sono assegnati da fuori ma calcolati
        K
        FF
        Cd0
        dCd0
        Re
        Cf
        e
        % Ipersostentatori
        yob                                     % Posizione della Sezione dell'ipersost. lungo y [m]
        cfoc                                    % Rapporto corda ipersost. corda profilo [-]
        cextoc                                  % Rapporto corda profilo con ipersost. esteso su corda senza ipersost. [-]
        flaptype
        auxvals                                 %variabile che contiene tutti i valori ausiliari nel calcolo delle grandezze
        deltaS                                  % Deflessioen degli Slats [deg]
        deltaF
        % Dati Equilibratore
        tauE                                    % Potenza di controllo dell'equilibratore
        eps0
        depsda
        % Global Positions
        xglob
        yglob
        zglob
    end
    methods
        function obj = ProfileClass(geomV,aeroV,M,HLflag,HLvals) %Costruttore
            %   geomV: array contenente le informazioni geometriche
            %   della sezione
            %   aeroV: array che contiene le informazioni aerodinamiche
            %       delal sezione organizzate per righe: una riga corrisponde
            %   a   d una condiziopne di volo.
            %   M: mach di volo al quale sono definite le proprietà
            %       aerodinamiche
            %   HLflag: variabile che assume 1: se il dato di HLvalòs si
            %       riferisce a flap; 2 se slat; 3 se piano orizzontale
            %   HLvals: vvalori di slat e/o flap oppure del piano di coda
            obj.K =NaN; obj.Cd0 =NaN; obj.Re =NaN; obj.Cf =NaN; obj.e = 1;
            obj.FF =NaN; obj.h = NaN;
            % Assegna le coordinate del profilo pulito
            obj.c = geomV(1);
            obj.eps = geomV(2);
            obj.tc = geomV(3);
            obj.LERc = geomV(4);
            obj.xtc = geomV(5);
            obj.xtrUp = geomV(6);
            obj.xtrLow = geomV(7);
            obj.dY = geomV(8);
            if length(M) == length(aeroV(:,1))
                obj.M = M;
                obj.a = aeroV(:,1);
                obj.cl0 = aeroV(:,2);
                obj.clstar = aeroV(:,3);
                obj.clmax = aeroV(:,4);
                obj.alphamax = aeroV(:,5);
                obj.alpha0l = aeroV(:,6);
                obj.alphastar = aeroV(:,7);
                obj.cmac = aeroV(:,8);
            else
                error("Numero di variabili aerodinamiche diverso dal numero di mach dati")
            end
            % Inizializza con NaN tutte le variabili non assegnate
            % dall'esterno
            obj.deltaF = NaN; obj.deltaS = NaN;
            obj.xglob = NaN; obj.yglob = NaN;
            obj.zglob = NaN; obj.yob = NaN;
            obj.flag = 'Panel Profile';
            obj.cfoc = NaN; obj.flaptype = NaN;
            obj.cextoc = NaN; obj.tauE = NaN;
            obj.eps0 = NaN; obj.depsda = NaN;
            obj.auxvals = NaN(10,1);
            if nargin == 5
                switch HLflag
                    case 'flaps'
                        obj.flag = 'Flapper Profile';
                        obj.yglob = HLvals(1);
                        obj.yob = HLvals(2);
                        obj.cfoc = HLvals(3);
                        switch HLvals(4)
                            case 1
                                obj.flaptype = 'fowler';
                            case 2
                                obj.flaptype = 'plain';
                            case 3
                                obj.flaptype = '3-slotted';
                            case 4
                                obj.flaptype = '2-slotted';
                        end
                    case 'slats'
                        obj.flag = 'Slat Profile';
                        obj.yglob = HLvals(1);
                        obj.yob = HLvals(2);
                        obj.cfoc = HLvals(3);
                        obj.cextoc = HLvals(4);
                    case 'elevator'
                        obj.flag = 'Horizontal Tail';
                        obj.tauE = HLvals(1);
                        obj.eps0 = HLvals(2);
                        obj.depsda = HLvals(3);
                end
            end
        end
        % function [X,Y,Z] = globalPos(obj,X,Y,Z)
        %     if nargin < 3
        %         Z = 0;
        %     end
        % end
        % function r = multiplyBy(obj,n)
        %     r = [obj.Value]*n;
        % end
        function out = HLauxVariables(obj,flag,var)
            if nargin == 3
                switch flag
                    case 'dCl02D'
                        obj.auxvals(1) = var;
                    case 'dClmax2D'
                        obj.auxvals(2) = var;
                    case 'a'
                        obj.auxvals(3) = var;
                    case 'alphaDeltaf'
                        obj.auxvals(4) = var;
                    case 'cbaroc'
                        obj.auxvals(5) = var;
                    case 'foo'
                        obj.auxvals(6) = var;
                        % Grandezze per Slats
                    case 'etaMax'
                        obj.auxvals(7) = var;
                    case 'etaDelta'
                        obj.auxvals(8) = var;
                    case 'CloDs'
                        obj.auxvals(9) = var;
                    case 'dClmax_slat'
                        obj.auxvals(10) = var;
                    otherwise
                        error('Campo inesistente')
                end
            else
                switch flag
                    case 'dCl02D'
                        out = obj.auxvals(1);
                    case 'dClmax2D'
                        out = obj.auxvals(2);
                    case 'a'
                        out = obj.auxvals(3);
                    case 'alphaDeltaf'
                        out = obj.auxvals(4);
                    case 'cbaroc'
                        out = obj.auxvals(5);
                    case 'foo'
                        out = obj.auxvals(6);
                        % Grandezze per Slats
                    case 'etaMax'
                        out = obj.auxvals(7);
                    case 'etaDelta'
                        out = obj.auxvals(8);
                    case 'CloDs'
                        out = obj.auxvals(9);
                    case 'dClmax_slat'
                        out = obj.auxvals(10);
                    otherwise
                        error('Campo inesistente')
                end
            end
        end
        function CL = CLvsAlphaCurve(obj,alpha)
        %CLvsAlphaCurve restituisce la curva di portanza del profilo al
        %               variare di alpha

            % CL = (alpha*obj.a + obj.cl0).*(alpha<obj.alphastar) + ...
            %     interp1([obj.alpha0l,obj.alphastar,obj.alphamax],[0,obj.clstar,obj.clmax],alpha,'spline',0).*(~(alpha<obj.alphastar));
            % Costruzione della Matrice dei coefficienti della cubica:
            % cfs = [a b c] a + b*x + c*x^2
            Mt = [1 obj.alphastar obj.alphastar^2 obj.alphastar^3; ...
                0 1 2*obj.alphastar 3*obj.alphastar^2; ...
                1 obj.alphamax obj.alphamax^2 obj.alphamax^3;...
                0 1 2*obj.alphamax 3*obj.alphamax^2];
            tn = [obj.clstar,obj.a,obj.clmax,0];
            cfs = Mt\tn(:);
            CL = (alpha*obj.a + obj.cl0).*(alpha<obj.alphastar) + ...
                ( cfs(1) + cfs(2)*alpha + cfs(3)*alpha.^2  + cfs(4)*alpha.^3 ).*(~(alpha<obj.alphastar));
            

        end
    end
end