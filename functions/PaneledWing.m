classdef PaneledWing
   properties
       m
       M
       geom_sect
       int_sect
   end
   methods
       function obj = PaneledWing(m,M,geom_vec,aero_vec,b)
           %m: number of sections
           %M: number of integration points
           %geom_vec: vector of geometric sections 
           %    c,eps,t/c,LER/c,x@max(t/c),xtr_Up,xtr_Low,dY
           %aero_vec: vector of aerodynamic data for profiles:
           %    Mach,cla,cl0,cl*,clmax,alphamax,alpha0l,alpha*,cm_ac
           % TODO: given c_vec, interpolate at y
           obj.m = m;
           obj.M = M;
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

       
   end
end