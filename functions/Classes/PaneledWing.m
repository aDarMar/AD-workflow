classdef PaneledWing < WingClass
   properties
       m        % number of spanwise sections
       m_red    % (m+1)/2
       M        % number of integration sections
       input_geom_sect
       geom_sect
       int_sect
       b
       sweep    % c/4 sweep [deg]
       % Global Profiles
       a_coeffs
       b_coeffs
       

   end
   methods
       function obj = PaneledWing(m,M,geom_vec,aero_vec,b,sweep,dihedral,iang,apexC,Mach)
           %m: number of sections
           %M: number of integration points
           %geom_vec: vector of half-wing geometric sections 
           %    2y/b,c,eps,t/c,LER/c,x@max(t/c),xtr_Up,xtr_Low,dY
           %aero_vec: vector of aerodynamic data for profiles:
           %    Mach,cla,cl0,cl*,clmax,alphamax,alpha0l,alpha*,cm_ac
           %    all angles are given in deg
           %    b: wing span
           %    sweep: leading edge sweep in degs
           % TODO: given c_vec, interpolate at y
           
           %% Superclass Construction 
           % Identification of Kink Section: checks at which station the
           % chord distribution is discontinous
           m_inp = length( geom_vec(:,1) );
           j = 1; kink_idx = nan( m_inp-1,1 );
           k = 1;
           % Calculates the slope of c - y for k = 1, k = 2
           % tg_old = c(2) - c(1) / y(2) - y(1)
           tg_old = ( geom_vec(k+1,2) - geom_vec(k,2) )/( geom_vec(k+1,1) - geom_vec(k,1) )/(0.5*b);
           for k = 2:m_inp-1
               tg_k = ( geom_vec(k+1,2) - geom_vec(k,2) )/( geom_vec(k+1,1) - geom_vec(k,1) )/(0.5*b);
               if abs( tg_k - tg_old ) > 1e-2
                   % if two consecutive tangents are discontinous then
                   % there is a kink at k station
                  kink_idx(j) = k;
                  j = j+1;
               end
               tg_old = tg_k;
           end
           kink_idx = kink_idx( ~isnan(kink_idx) ); % saves index of kink positions
           n_pan    = length( kink_idx ) + 1; % number of panels = n_kink + 1
           kink_idx = [1;kink_idx;m_inp]; % indexes of sections of interest [ root, kink_1,kink_2,..,tip]
           % Semi-span evaluation
           bs = nan(n_pan,1);
           for k = 1:n_pan
               bs(k) = ( geom_vec( kink_idx(k+1),1 ) - ...
                   geom_vec( kink_idx(k),1 ) )*b/2;
           end
           % WARNING: only one dihedral and sweep for the entire wing
           obj@WingClass( bs,ones(n_pan)*sweep,ones(n_pan)*dihedral,iang,apexC,Mach,...
               geom_vec(kink_idx,2:end),aero_vec(kink_idx,2:end) ); % call superclass constructor
           
           %% Sections Definition

           if mod( m,2 ) == 0
               warning( 'Inserted even numer of points, added one point to make it odd' );
               m = m+1;
           end
           % Wing Sections
           % x -- o -- o -- |o| -- o -- o -- x
           % |    m       (m+1)/2  2    1    |
           % |                               |
           % |<------------- b ------------->|
           obj.m     = m; obj.m_red = ( obj.m + 1 ) *0.5; % Half points
           obj.M     = M;
           obj.b     = b;
           obj.sweep = obj.sweepChange( sweep,0,0.25,obj.AR,obj.TR );
           phi       = obj.phi_funct(obj.m);
           % Stores the input wing sections
           obj.input_geom_sect  = ProfileClass.empty; % Costruisce un array di oggetti panels
           for i=1:m_inp
               obj.input_geom_sect(i)       = ProfileClass( geom_vec(i,2:end),...
                   aero_vec(i,2:end),aero_vec(i,1) );
               obj.input_geom_sect(i).yglob = geom_vec(i,1)*obj.b/2;
           end
           % WARNING: in this section angles are defined in radiants
           % because og the calculations involved i nWeissinger method,
           % while in the previous section were given in degrees.
           
           % Defines the wing sections employed in calculations
           [geom_prep,aero_prep] = obj.interpSects(geom_vec,aero_vec,phi); % angles ae converted in radiants
           obj.geom_sect  = panelSection.empty; % Costruisce un array di oggetti panels
           for i=1:obj.m
               obj.geom_sect(i) = panelSection( geom_prep(i,:),aero_prep(i,2:end),aero_prep(i,1),phi(i),b );
           end
           % Defines the integration points sections
           obj.int_sect  = panelSection.empty; % Costruisce un array di oggetti panels
           temp = nan(1,8);
           phi   = obj.phi_funct(obj.M);
           for i=1:obj.m
               obj.int_sect(i) = panelSection( temp,temp,temp(1),phi(i),b );
           end
           % Plantform wing parameters
           %            obj.TR    = obj.input_geom_sect(1).c/obj.input_geom_sect( end ).c;
           %            obj.S     = obj.areacalc;
           %            obj.AR    = obj.b^2/obj.S;
           %            obj.sweepLE = obj.sweepChange(obj.sweep,0.25,0,obj.AR,obj.TR);
           %% Weissinger
           % Calculates Wing Loading
           [obj.a_coeffs,obj.b_coeffs] = obj.aeroDef; % Defines influence ad geometric coefficients
           [ Gb,Ga,CLa,a0L,CM0_basic ] = obj.loadsEval;
           for n = 1:obj.m_red-1
               % Symmetric Loadings
               obj.geom_sect(n).Gb = Gb(n); obj.geom_sect(obj.m+1-n).Gb = Gb(n);
               obj.geom_sect(n).Ga = Ga(n); obj.geom_sect(obj.m+1-n).Ga = Ga(n);
           end
           n = obj.m_red;
           obj.geom_sect(n).Gb = Gb(n);
           obj.geom_sect(n).Ga = Ga(n);
       end
       
       function [geom_prep,aero_prep] = interpSects( obj,geom_vec,aero_vec,phi )
           %interpSects: function that interpolates the characteristic from
           %input sections at sections used by the Multhopp integration
           %formula
           %OUTPUT
           %    aero_prep: aerodynamic data of the sections given in
           %    [rad]and [1/rad]. <- BEWARE!
           geom_prep = nan( obj.m,8 ); aero_prep = nan( obj.m,9 );
           for n = 1:obj.m
               for l = 1:8
                  geom_prep(n,l) = interp1( geom_vec(:,1),geom_vec(:,l+1),cos( phi(n) ) );
                  aero_prep(n,l) = interp1( geom_vec(:,1),aero_vec(:,l),cos( phi(n) ) );
               end
           end
           n = 2; geom_prep(:,n) = geom_prep(:,n)*pi/180;   % eps in [rad]
           aero_prep(:,2) = aero_prep(:,2)*180/pi;          % cla in [1/rad]
           aero_prep(:,6:8) = aero_prep(:,6:8)*pi/180;      % alphas in [rad]
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
       
       function Sw = areacalc(obj)
           %areacalc: calculates wing area from wing slices
           m_inp = length( obj.input_geom_sect );
           Sw = 0;
          for i = 2:m_inp
              Sw =  Sw + (obj.input_geom_sect(i-1).c+obj.input_geom_sect(i).c)*(obj.input_geom_sect(i-1).yglob-obj.input_geom_sect(i).yglob);
          end
       end
       
       function phi_n = phi_funct( ~,m )
           %phi_funct: function that calculates spanwise sections according to
           %Multhopp quadrature method
           %   m: number of points
           phi_n = nan(m,1);
           for n = 1:m
               phi_n(n) = n*pi/(m+1);
           end
       end
       
       function eps = eps_calc(obj,nu_r)
           %eps_calc: calculates geometric and aerodynamic twist for all
           %the sections in relation to the reference (root) section
           %INPUT
           %    nu_r: index of reference station
           %OUTPUT
           %    eps: vector of combined twist, given from (looking towards 
           %        the wing's apex) right tip to root to left tip
          eps = nan( obj.m_red,1 ); 
          for n = 1:obj.m_red
              eps(n) = obj.geom_sect(n).eps - obj.geom_sect(nu_r).eps; %+ % FINIREEEEEEEEE
          end
          
       end
       %% Aerodynamic Estimation
       function [a_cfs,b_cfs] = aeroDef(obj)
           a_cfs        = nan(obj.m_red);       %Matrix of symmetric influence coefficients
           b_cfs        = nan(obj.m_red);
           for nu = 1:obj.m_red
               % Builds the matrix of coefficients row by row
               [a_cfs(nu,:),b_cfs(nu,:)] = obj.aerosymmbuilder( nu );
           end
       end
       
       %% Influence Coefficients Definition
       function [a_nu,B] = aerosymmbuilder( obj,nu )
           %aerosymmbuilder: function that builds the equation for circulation in case of
           %symmetric load. In other words, it builds the nu-th row of the influence
           %coefficients matrix.
           %INPUT
           %   nu: nu-th row of equation. In physical terms it represents the index of
           %       the nu-th control point in which we are imposing the no flow
           %       through condition. NU IS A SCALAR
           %OUTPUT
           %    a_nu: matrix of influence coefficients: 
           %            [a_nu]{Gs} = {alpha_nu}
           %    B: matrix of geometric coefficients
           
           a_nu    = nan( 1,obj.m_red );
           n_idxs  = 1:obj.m_red-1;
           n_idxs  = n_idxs( n_idxs~=nu );
           B       = nan( 1,obj.m_red);
           % For-cycle from 1 to (m+1)/2 - 1 excluding nu
           for n = n_idxs
               % Eq. (A37) case n =/= nu
               B(n)    = obj.littlebfun( nu,n ) + obj.littlebfun( nu,obj.m+1-n );
               a_nu(n) = -2*B(n) + obj.b/( obj.geom_sect(nu).c*obj.geom_sect(nu).cor_2pi )*obj.gbarfun( nu,n );
           end
           % Adding (m+1)/2 point: this step is made outside the for cycle because the
           % functions gbar and B assume  different values for n = (m+1)/2
           B(obj.m_red)    = obj.littlebfun( nu,obj.m_red );
           a_nu(obj.m_red) = -2*B(obj.m_red) + obj.b/( obj.geom_sect(nu).c*obj.geom_sect(nu).cor_2pi )*obj.gbarfun_special( nu,obj.m_red );
           % Adding nu point
           B(nu)       = (obj.m+1)/( 4*sin(obj.geom_sect(nu).phi) ); % b(nu,nu)
           a_nu(nu)    = 2*B(nu);    %2*b_nu,nu expression
           if nu == obj.m_red
               % case nu = n = (m+1)/2
               a_nu(nu) = a_nu(nu) + obj.b/( obj.geom_sect(nu).c*obj.geom_sect(nu).cor_2pi )*obj.gbarfun_special( nu,nu );
           else
               % case nu = n =/= (m+1)/2
               a_nu(nu) = a_nu(nu) + obj.b/( obj.geom_sect(nu).c*obj.geom_sect(nu).cor_2pi )*obj.gbarfun( nu,nu );
           end
           
           
       end

       function b_ij = littlebfun( obj,nu,n )
           %littlebfun: function that evaluates the coefficient b(nu,n)
           %   nu: index of the control point
           %   n: index of the other point
           b_ij = sin( obj.geom_sect(n).phi )/( ( obj.geom_sect(n).eta-obj.geom_sect(nu).eta )^2 )*( (1-(-1)^(n-nu))/(2*(obj.m+1)) );
       end
       
       function gval = gbarfun( obj,nu,n )
           %gbarfun: function that evaluates g_bar in case the second index is
           %different from (m+1)/2.
           % ------------------ BE CAREFUL! ------------------
           %In case the value of g for that index is required,
           %refer to the function gbarfun_special
           %--------------------------------------------------
           %    nu: index of the control point
           %    n: index of the other point
           % Case mu = 0 and n =/= (m+1)/2
           gval = obj.ffun( obj.geom_sect(n).phi,0 )*...
               obj.L_funsymm( nu,1 );
           Mred = (obj.M-1)/2;
           for mu = 1:Mred
               % Case mu =/= 0 and n =/= (m+1)/2
               gval = gval + 2*obj.ffun( obj.geom_sect(n).phi,obj.int_sect(mu).phi )...
                   *obj.L_funsymm( nu,obj.int_sect(mu).eta );
           end
           gval = gval*(-1)/(2*(obj.M+1));
       end
       
       function gval = gbarfun_special( obj,nu,n )
           %gbarfun: function that evaluates g_bar in case the second index IS
           %EQUAL TO (m+1)/2.
           % ------------------ BE CAREFUL! ------------------
           %In case the value of g for the other indexes is required,
           %refer to the function gbarfun
           %--------------------------------------------------
           %    nu: index of the control point
           %    n: index of the other point
           
           % Case mu = 0 and n == (m+1)/2
           gval = 0.5*obj.ffun( obj.geom_sect(n).phi,0 )*...
               obj.L_funsymm( nu,1 );
           Mred = (obj.M-1)/2;
           for mu = 1:Mred
               % Case mu =/= 0 and n == (m+1)/2
               gval = gval + obj.ffun( obj.geom_sect(n).phi,obj.int_sect(mu).phi )...
                   *obj.L_funsymm( nu,obj.int_sect(mu).eta );
           end
           gval = gval*(-1)/(2*(obj.M+1));
       end
       
       function fnm = ffun( obj,phi_n,phi_mu )
           %    phi_n: geometrical section index
           %    phi_mu: integration point index
           fnm = 0;
           idx = 1:2:obj.m;
           for mu1 = idx
               fnm = fnm + mu1*sin( mu1*phi_n )*cos( mu1*phi_mu );
           end
           fnm = fnm*2/(obj.m+1);
       end
       
       function L = L_funsymm( obj,nu,b_eta ) % CONTROLLA LE CHIAMATEEEEEEE
           %L_funsymm: Weissner's influence function for a symmetric load 
           %distribution. The reason the function requires the index for the 
           %control point and the actual position for the integration point 
           %is that the second index can be 0, while the first cant.
           %    nu: index of control point nu
           %    b_eta: non-dimensional position of intehration point mu
           boc = obj.b/( obj.geom_sect(nu).c*obj.geom_sect(nu).cor_2pi ); tS4 = tan(obj.sweep*pi/180)/obj.geom_sect(nu).beta;
           eta = obj.geom_sect(nu).eta;
           if abs(eta-b_eta) < 1e-4
               %If eta = b_eta it can be shown that L 
               L = tS4;
           else
               L = 1/( boc*(eta-b_eta) )*( sqrt( (1+boc*(eta-b_eta)*tS4)^2 + (boc*(eta-b_eta))^2 ) - 1 );
           end
           L = L - 1/( boc*(eta+b_eta) )*( sqrt( (1+boc*(eta-b_eta)*tS4)^2 + (boc*(eta+b_eta))^2 )/( 1+2*boc*eta*tS4 ) - 1 ) -...
                   ( 2*tS4*sqrt( (1+boc*eta*tS4)^2+(boc*eta)^2 ) )/(1+2*boc*eta*tS4);
       end
       
       %% Loads Calculation
       function [ Gb,Ga,CLa,a0L,CM0_basic ] = loadsEval(obj)
           [Gb, a0L] = obj.basicLoad;              % Basic load distribution
           %CDi_basic = obj.induced_drag( Gb );     % Induced drag due to Basic Loads
           CM0_basic = obj.pitchCoeff_basic( Gb ); % Pitching Moment Coeff. due to Basic load
           
           [Ga,CLa] = obj.additionalLoad;           % Gna/alpha, CLa [1/rad]
           %CDi_add  = obj.induced_drag( Ga );       % Induced drag due to Basic Loads
           
           %Gtot = Ga + Gb;
           %CDi_tot = obj.induced_drag( Gtot );
       end
       
       function CDi = induced_drag( obj,Gs )
           %induced_Drag: function that evaluates partial induced drag 
           %coefficient (only basic or additional)
           %INPUT
           %    Gs: circulation along stations from root to centerline
           CDi = 0;
           m_rd = obj.m_red - 1; idx_b = 1:obj.m_red; 
           for nu = 1:m_rd
               temp = 0;
               idx = idx_b( idx_b ~= nu );
               for n = idx
                   % Summing the B(nu,n)G(n)
                   temp = temp + obj.b_coeffs(nu,n)*Gs(n);
               end
               CDi = CDi + 2*( obj.b_coeffs(nu,nu)*Gs(nu) - temp )*sin( obj.geom_sect(nu).phi )*Gs(nu);
           end
           nu = idx_b(end); temp = 0;
           idx = idx_b( idx_b ~= nu );
           for n = idx
               temp = temp + obj.b_coeffs(nu,n)*Gs(n);
           end
           CDi = CDi + ( obj.b_coeffs(nu,nu)*Gs(nu) - temp )*sin( obj.geom_sect(nu).phi )*Gs(nu);
           CDi = CDi*pi*obj.AR/(obj.m+1);
       end
       
       function CL = liftCoeff( obj,G )
           CL = 0; m_rd = obj.m_red-1;
           for n = 1:m_rd
              CL = CL + 2*G(n)*sin( obj.geom_sect(n).phi );
           end
           CL = CL + G(obj.m_red); CL = CL*obj.AR*pi/8;
       end
       % Basic Load
       function [Gs,a0L_ref] = basicLoad( obj,nu_r )
           %basicLoad: function that evaluates the circulation due to the
           %basic load
           %INPUT
           %    nu_r: station chosen as reference
           if nargin < 2
               nu_r   = obj.m_red; % reference station is the root chord
           end
           eps = obj.eps_calc( nu_r );
           % Reduced Matrix for Basic Loading
           % i: i-th control point (nu in reference)
           % j: j-th element
           % r: index of reference point. IT IS ASSUMED THIS IS THE ROOT, otherwise
           % the equation below is no longer valid
           % A(i,j) = a(i,j)-a(r,j)-( a(i,r)-a(r,r) )*2*sin( phi(j) )
           
           nu_idx = 1:obj.m_red; nu_idx_red = nu_idx( nu_idx ~= nu_r);
           A      = nan( obj.m_red - 1 );
           for nu = nu_idx_red
               % Twist in radiants
               eps(nu) = ( obj.geom_sect(nu).eps + obj.geom_sect(nu).alpha0l - obj.geom_sect(nu_r).alpha0l )*pi/180;
               for n = nu_idx_red
                   A(nu,n) = obj.a_coeffs(nu,n)-obj.a_coeffs(nu_r,n)-...
                       ( obj.a_coeffs(nu,nu_r)-obj.a_coeffs(nu_r,nu_r) )*2*sin( obj.geom_sect(n).phi );
               end
           end
           Gs      = A\eps(nu_idx_red);
           if nu_r > 1
               Gs  = [Gs(1:nu_r-1);0;Gs(nu_r:end)];
           else
               Gs  = [0;Gs];
           end
           a0L_ref = 0;
           for n = nu_idx_red
               % Equation for total CL = 0
               Gs(nu_r) = Gs(nu_r) - Gs(n)*2*sin( obj.geom_sect(n).phi );
               % Equation for circulation at Reference Station
               a0L_ref  = a0L_ref + ( obj.a_coeffs(nu_r,n) - obj.a_coeffs(nu_r,nu_r)...
                   *2*sin( obj.geom_sect(n).phi ) )*Gs(n);
           end
       end
       
       function Cm0 = pitchCoeff_basic( obj,Gs )
           %basicLoad_coeffs: function that evaluates the pitching moment
           %due to the basic load using eq (10) of NACA TR921. For details
           %about the formulas see the doc.
           %INPUT
           %    Gs: vector of spanwise circulations, from b/2 to 0
           A = zeros( 1,obj.m_red );
           oddx = 1:2:obj.m;	 %odd indices for mu_1 because of the integral sin(mu*phi)sin(phi)
           j = 1;
           for n = 1:obj.m_red
              for mu = oddx
                  A(j) = A(j) + sin(mu*pi/2)*sin( mu*obj.geom_sect(n).phi )*2/(4-mu^2);
              end
              j = j + 1;
           end
           A(1:obj.m_red-1) = 2*A(1:obj.m_red-1); % Coefficients for stations 1 to the first before ceter span are doubled
           A = A/(obj.m+1);
           %cm0_basic = A*obj
           Cm0 = -obj.AR*obj.b*tan(obj.sweep*pi/180)/obj.meanprofile.c*A*Gs; %FINIREEEEEE
       end
       % Additional Load
       function [Ga,CLa] = additionalLoad( obj )
           RHS = ones( (obj.m+1)*0.5,1 );
           Ga = obj.a_coeffs\RHS;
           CLa = obj.liftCoeff( Ga );
       end
       
       % Plot
       function wing_circ(obj,alpha)
           alpha = alpha*pi/180;
            fig    = figure("Name",'Wing Span Load'); 
            ax_fig = axes('Parent',fig);hold( ax_fig,'on' );
            Ga_v = nan(obj.m_red,1); Gb_v = Ga_v; eta = Ga_v; Gtot = Ga_v;
            for i = 1:obj.m_red
                Ga_v(i) = obj.geom_sect(i).Ga;
                Gb_v(i) = obj.geom_sect(i).Gb;
                eta(i)  = obj.geom_sect(i).eta;
                Gtot(i) = Ga_v(i)*alpha + Gb_v(i);
            end
            Gtot = [0;Gtot]; eta = [1;eta];
            Gb_v = [0;Gb_v]; Ga_v = [0;Ga_v];
            lin(1) = plot( ax_fig,eta,Gtot );
            lin(2) = plot( ax_fig,eta,Gb_v );
            lin(3) = plot( ax_fig,eta,Ga_v*alpha );
       end
       % NON-Weissinger
       % CM and Alpha0L
   end
   
end