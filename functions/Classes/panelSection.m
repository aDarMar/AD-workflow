classdef panelSection < ProfileClass
    properties
        % Nondimensional section coordinates:
        % y = b/2 * eta = b/2 * cos( phi )
        eta
        phi
        % Nondimensional loadings:
        % G = c(y)écl(y)/(2*b)
        Gb   % Basic
        Ga   % Additional
        G    % Total
        % Corrections
        cor_2pi     % Ration between actual and ideal lift slope (radiants)
        beta        % Prandtl-Glauert Compressibility factor
   end
    methods
        function obj = panelSection(geomV,aeroV,M,phi,b)
           obj@ProfileClass(geomV,aeroV,M) % call superclass constructor
           obj.phi   = phi;
           obj.eta   = cos(obj.phi);
           obj.yglob = obj.eta*b*0.5;
           if obj.M ~= 0 && ~isnan( obj.M )
               if obj.M > 1
                   obj.beta = sqrt( obj.M^2-1 );
               else
                  obj.beta = sqrt( 1-obj.M^2 ); 
               end
           else
               obj.beta = 1;
           end
           if obj.a ~= 0 && ~isnan( obj.a )
               % There is a double correction: one for the real lift slope
               % that is cl_a/(2*pi/b) and the other is for the
               % compressibility 1/b and are combined in such a way that b
               % cancels out
               obj.cor_2pi = obj.a/(2*pi);
           else
               obj.cor_2pi = 1;
           end
        end
    end
    
    
    
end