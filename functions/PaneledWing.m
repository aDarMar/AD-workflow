classdef PaneledWing
   properties
       m
       M
       geom_sect
       int_sect
       b
       sweep    % c/4 sweep [deg]
       sweepLE  % LE sweep
       TR
       S
       AR
   end
   methods
       function obj = PaneledWing(m,M,geom_vec,aero_vec,b,sweep)
           %m: number of sections
           %M: number of integration points
           %geom_vec: vector of geometric sections 
           %    c,eps,t/c,LER/c,x@max(t/c),xtr_Up,xtr_Low,dY
           %aero_vec: vector of aerodynamic data for profiles:
           %    Mach,cla,cl0,cl*,clmax,alphamax,alpha0l,alpha*,cm_ac
           %    b: wing span
           %    sweep: leading edge sweep in degs
           % TODO: given c_vec, interpolate at y
           if mod( m,2 ) == 0
               warning( 'Inserted even numer of points, added one point to make it odd' );
               m = m+1;
               % FINIREEEEEEEEEE: aggiungere sezione mancante
           end
           obj.m = m;
           obj.M = M;
           obj.b = b;
           obj.sweep = sweep;
           phi   = obj.phi_funct(obj.m);
           % Defines the wing sections
           obj.geom_sect  = ProfilelLass.empty; % Costruisce un array di oggetti panels
           for i=1:obj.m
               obj.geom_sect(i) = panelSection( geom_vec(i,:),aero_vec(i,2:end),aero_vec(i,1),phi(i),b );
           end
           % Defines the integration points sections
           obj.int_sect  = ProfilelLass.empty; % Costruisce un array di oggetti panels
           temp = nan(1,8);
           phi   = obj.phi_funct(obj.M);
           for i=1:obj.m
               obj.int_sect(i) = panelSection( temp,temp,temp(1),phi(i),b );
           end
           % Plantform wing parameters
           obj.TR    = obj.geom_sect(1).c/obj.geom_sect( (obj.m+1)*0.5 ).c;
           obj.S     = obj.areacalc;
           obj.AR    = obj.b^2/obj.S;
           obj.sweep = sweepChange(obj.sweepLE,0,0.25,obj.AR,obj.TR);
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
           Sw = 0;
          for i = 2:(obj.m+1)*0.5 
              Sw =  Sw + (obj.geom_sect(i-1).c+obj.geom_sect(i).c)*0.5*obj.b*(obj.geom_sect(i-1).eta-obj.geom_sect(i).eta);
          end
       end
       
       function phi_n = phi_funct( ~,m )
           %phi_funct: function that calculates spanwise sections according to
           %Multhopp quadrature method
           %   m: number of points
           phi_n = nan(m,1);
           for i = 1:m
               phi_n(i) = n*pi/(m+1);
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
               obj.L_funsymm( obj.geom_sect(nu).eta,1 );
           Mred = (obj.M-1)/2;
           for mu = 1:Mred
               % Case mu =/= 0 and n =/= (m+1)/2
               gval = gval + 2*obj.ffun( obj.geom_sect(n).phi,obj.int_sect(mu).phi )...
                   *obj.L_funsymm( obj.geom_sect(nu).eta,obj.int_sect(mu).eta );
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
               obj.L_funsymm( obj.geom_sect(nu).eta,1 );
           Mred = (obj.M-1)/2;
           for mu = 1:Mred
               % Case mu =/= 0 and n == (m+1)/2
               gval = gval + obj.ffun( obj.geom_sect(n).phi,obj.int_sect(mu).phi )...
                   *obj.L_funsymm( obj.geom_sect(nu).eta,obj.int_sect(mu).eta );
           end
           gval = gval*(-1)/(2*(obj.M+1));
       end
       
       function fnm = ffun( obj,phi_n,phi_mu )
           %    phi_n: geometrical section index
           %    phi_mu: integration point index
           for mu1 = 1:obj.m
               fnm = fnm + mu1*sin( mu1*phi_n )*cos( mu1*phi_mu );
           end
           fnm = fnm*2/(obj.m+1);
       end
       
       function L = L_funsymm( obj,eta,b_eta ) % CONTROLLA LE CHIAMATEEEEEEE
           %L_funsymm: Weissner influence function for a symmetric load
           %distribution.
           %    eta: non-dimensional position of control point nu
           %    b_eta: non-dimensional position of intehration point mu
           boc = obj.b/obj.c; tS4 = tan(bj.sweep*pi/180);
           L = 1/( boc*(eta-b_eta) )*( sqrt( (1+boc*(eta-b_eta)*tS4)^2 + (boc*(eta-b_eta)^2) ) - 1 ) -...
               1/( boc*(eta-b_eta) )*( sqrt( (1+boc*(eta-b_eta)*tS4)^2 + (boc*(eta-b_eta)^2) )/( 1+2*boc*eta*tS4 ) - 1 ) -...
               ( 2*tS4*sqrt( (1+boc*eta*tS4)^2+(boc*eta)^2 ) )/(1+2*boc*eta*tS4);
       end
       
   end
end