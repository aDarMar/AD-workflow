classdef PanelClass < handle%< ProfileClass
    properties
        % Panel Geometry
        b {mustBeNumeric, mustBePositive}             %apertura [m]
        sweep {mustBeNumeric}                         %freccia [deg]
        dihedral
        S {mustBeNumeric, mustBePositive}             %Superficie (m^2)
        mac {mustBeNumeric, mustBePositive}           %Corda media aerodinamica(m)
        root
        tip
        TR {mustBeNumeric, mustBePositive}            %Taper Ratio [-]
        xmac
        ymac                                          %Coordinate della corda media aerodinamica nel rif. pannello
        zmac
        % Panel High-Lift
        auxvals
        deltaCoeffs
        HLflap
        HLslat
    end
    methods
        function obj = PanelClass(b,Sweep,Dihedral,M,rootGeom,rootAero, ...
                tipGeom, tipAero, HLflag,HLvals) %Costruttore
            if nargin == 1
                %riceve l'indirizzo del file xml
            else


                if nargin < 9
                    % High Lift Devices not present
                    obj.root = ProfileClass(rootGeom,rootAero,M);
                    obj.tip  = ProfileClass(tipGeom,tipAero,M);
                else
                    % High lift devices present
                    obj.root = ProfileClass(rootGeom,rootAero,M,HLflag,HLvals(1,:));
                    obj.tip  = ProfileClass(tipGeom,tipAero,M,HLflag,HLvals(2,:));
                end
            end
            obj.b        = b;
            obj.sweep    = Sweep;
            obj.dihedral = Dihedral;
            obj.TR       = obj.tip.c/obj.root.c;
            obj.S        = areaFun(obj);
            obj.mac      = macCalc(obj);
            [obj.xmac,obj.ymac,obj.zmac] = macCoordinates(obj);
            
            %auxvals = [dCl0,dClmax,a_medio,alphaDelta_medio,kb,kc
            obj.auxvals = NaN;
            obj.deltaCoeffs = NaN;
            % HLflap = [deltaCL0, deltaCLmax,
            obj.HLflap = NaN(7,1);
        end
        function [xroot,yroot,zroot] = globalCoords(obj,Xg,Yg,Zg)
            %globalCoords: function that evaluates the global coordinates
            %of the profile sections
            %INPUT
            %   Xg,Yg,Zg: global coordinate (ie in the costructive
            %       reference frame) of the inner LE section
            %OUTPUT
            %   xroot,ysec,zsec: global coordinates (ie in the costructive
            %       reference frame) of the section
            %xroot = Xg;
            xroot  = Xg + obj.b*tan(obj.sweep*pi/180);
            %yroot = Yg;
            yroot  = Yg + obj.b;
            %zroot = Zg;
            zroot  = Zg + obj.b*tan(obj.dihedral*pi/180);

        end
        function S = areaFun(obj)
            S = (obj.root.c + obj.tip.c )*obj.b*0.5;
        end
        function mac = macCalc(obj)
            mac =  2/3 * obj.root.c * (1+obj.TR+obj.TR^2)/(1 + obj.TR);
        end
        function [Xle , Yle, Zle] = macCoordinates(obj)
            Yle = obj.b/6 * (1 + 2*obj.TR)/(1 + obj.TR);
            Xle = Yle * tan(obj.sweep * pi/180);
            Zle = Yle * tan(obj.dihedral * pi/180);
        end

        function [temp,y] = panelInterp(obj,y,flg)
            % panelInterp: interpola le grandezze geometriche dei profili
            % partendo dai valori di corda e radice.
            %y: coordinata y globale
            % sectsFlaps = [y, cf/c, flaptype]
            % sectsSlats = [y, cs/c]
            % Deve essere UNA SOLA RIGA del vettore delle coordinate
            %^ f
            %|
            %|                               x
            %|         x                     |
            %|         |                     |
            %|         |                     |
            %-----------------------------------------------> y
            %      root.yglobal       root.yglobal+panel.b
            %          |-------panel.b-------|
            y = y(:); lgc = NaN(length(y),1);
            for j=1:length(y)
                lgc(j) = y(j)-obj.root.yglob>0 && y(j)-obj.tip.yglob<0;
            end
            y = y(lgc==1); %prende solo i valori di y compresi nel pannello

            if nargin == 3
                temp = NaN(length(y),1);
            switch flg
                case 1
                    temp(:,1) = ... % corda
                        (obj.tip.c - obj.root.c)/obj.b * (y - obj.root.yglob) + obj.root.c; % Interpolazione corda
                case 2
                    temp(:,1) = ... % epsilon
                        (obj.tip.eps - obj.root.eps)/obj.b * (y - obj.root.yglob) + obj.root.eps;
                case 3
                    temp(:,1) = ...
                        (obj.tip.tc - obj.root.tc)/obj.b * (y- obj.root.yglob) + obj.root.tc;
                case 4
                    temp(:,1) = ...
                        (obj.tip.LERc - obj.root.LERc)/obj.b * (y - obj.root.yglob) + obj.root.LERc;
                case 5
                    temp(:,1) = ...
                        (obj.tip.xtc - obj.root.xtc)/obj.b * (y - obj.root.yglob) + obj.root.xtc;
                case 6
                    temp = NaN(length(y)*length(obj.tip.a),1);
                    temp(:,1) = ... %a
                        (obj.tip.a - obj.root.a)./obj.b .* (y - obj.root.yglob) + obj.root.a;
                case 7
                    temp = NaN(length(y)*length(obj.tip.a),1);
                    temp(:,1) = ... %cl0
                        (obj.tip.cl0 - obj.root.cl0)./obj.b .* (y - obj.root.yglob) + obj.root.cl0;
                case 8
                    temp = NaN(length(y)*length(obj.tip.a),1);
                    nl = length(y);
                    nM = length(obj.root.a);
                    for i=1:nl
                        temp((i-1)*nM + 1:i*nM,1) = ... %a
                            (obj.tip.a - obj.root.a)./obj.b .* (y(i) - obj.root.yglob) + obj.root.a;
                        temp((i-1)*nM + 1:i*nM,2)  = ... %cl0
                            (obj.tip.cl0 - obj.root.cl0)./obj.b .* (y(i) - obj.root.yglob) + obj.root.cl0;
                        temp((i-1)*nM + 1:i*nM,3) = ... %cl0
                            (obj.tip.clstar - obj.root.clstar)./obj.b * (y(i) - obj.root.yglob) + obj.root.clstar;
                        temp((i-1)*nM+ 1:i*nM,4) = ... %clmax
                            (obj.tip.clmax - obj.root.clmax)./obj.b * (y(i) - obj.root.yglob) + obj.root.clmax;
                        temp((i-1)*nM + 1:i*nM,5) = ... %alphamax
                            (obj.tip.alphamax - obj.root.alphamax)./obj.b * (y(i) - obj.root.yglob) + obj.root.alphamax;
                        temp((i-1)*nM + 1:i*nM,6) = ... %alpha0l
                            (obj.tip.alpha0l - obj.root.alpha0l)./obj.b * (y(i) - obj.root.yglob) + obj.root.alpha0l;
                        temp((i-1)*nM + 1:i*nM,7) = ... %alphastar
                            (obj.tip.alphastar - obj.root.alphastar)./obj.b * (y(i) - obj.root.yglob) + obj.root.alphastar;
                        temp((i-1)*nM+ 1:i*nM,8) = ... %cmac
                            (obj.tip.cmac - obj.root.cmac)./obj.b * (y(i) - obj.root.yglob) + obj.root.cmac;
                    end

            end

            else
                temp = NaN(length(y),8);
                temp(:,1) = ... % corda
                    (obj.tip.c - obj.root.c)/obj.b * (y - obj.root.yglob) + obj.root.c; % Interpolazione corda
                temp(:,2) = ... % epsilon
                    (obj.tip.eps - obj.root.eps)/obj.b * (y - obj.root.yglob) + obj.root.eps;
                temp(:,3) = ...
                    (obj.tip.tc - obj.root.tc)/obj.b * (y- obj.root.yglob) + obj.root.tc;
                temp(:,4) = ...
                    (obj.tip.LERc - obj.root.LERc)/obj.b * (y - obj.root.yglob) + obj.root.LERc;
                temp(:,5) = ...
                    (obj.tip.xtc - obj.root.xtc)/obj.b * (y - obj.root.yglob) + obj.root.xtc;
                temp(:,6) = ...
                    (obj.tip.xtrUp - obj.root.xtrUp)/obj.b * (y - obj.root.yglob) + obj.root.xtrUp;
                temp(:,7) = ...
                    (obj.tip.xtrLow - obj.root.xtrLow)/obj.b * (y - obj.root.yglob) + obj.root.xtrLow;
                temp(:,8) = ...
                    (obj.tip.dY - obj.root.dY)/obj.b * (y - obj.root.yglob) + obj.root.dY;
            end

        end
    
        % function out = coeffsShiftFun(obj,flag,valAss)
        %     switch flag
        %         case 'dCLmax'
        %             obj.HLflap(2) = valAss;
        %     end
        % end
        function out = coeffsShiftFun(obj,flag,var)
            if nargin == 3
                switch flag
                    case 'dCL0'
                        obj.deltaCoeffs(1) = var;
                    case 'dCLmax'
                        obj.deltaCoeffs(2) = var;
                    case 'dCLmax_slat'
                        obj.deltaCoeffs(3) = var;

                end
            else
                switch flag
                    case 'dCL0'
                        out = obj.deltaCoeffs(1);
                    case 'dCLmax'
                        out = obj.deltaCoeffs(2);
                        % Grandezze Slat
                    case 'dCLmax_slat'
                        out = obj.deltaCoeffs(3);
                end
            end
        end

        function out = HLauxVariables(obj,flag,var)
            if nargin == 3
                switch flag
                    case 'dCl02D_mean'
                        obj.auxvals(1) = var;
                    case 'dCLmax_mean'
                        obj.auxvals(2) = var;
                    case 'a_mean'
                        obj.auxvals(3) = var;
                    case 'alphaDeltaf'
                        obj.auxvals(4) = var;
                    case 'Kb'
                        obj.auxvals(5) = var;
                    case 'Kc'
                        obj.auxvals(6) = var;
                    case 'cbaroc_mean'
                        obj.auxvals(7) = var;
                        % Grandezze Slat
                    case 'dCmax2D_slat'
                        obj.auxvals(8) = var;
                    otherwise
                        error('Campo non esistente');
                end
            else
                switch flag
                    case 'dCl02D_mean'
                        out = obj.auxvals(1);
                    case 'dCLmax_mean'
                         out = obj.auxvals(2);
                    case 'a_mean'
                        out = obj.auxvals(3);
                    case 'alphaDeltaf'
                        out = obj.auxvals(4);
                    case 'Kb'
                        out = obj.auxvals(5);
                    case 'Kc'
                        out = obj.auxvals(6);
                    case 'cbaroc_mean'
                        out = obj.auxvals(7);
                        % Grandezze Slat
                    case 'dCmax2D_slat'
                        out = obj.auxvals(8);
                    otherwise
                        error('Campo non esistente');
                end
            end
        end


    end

end
