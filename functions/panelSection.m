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
        cor_2pi
   end
    methods
        function obj = panelSection(geomV,aeroV,M,phi,b)
           obj@ProfileClass(geomV,aeroV,M) % call superclass constructor
           obj.phi   = phi;
           obj.eta   = cos(obj.phi);
           obj.yglob = obj.eta*b*0.5;
           if obj.a ~= 0 && ~isnan( obj.a )
               obj.cor_2pi = obj.a/(2*pi);
           else
               obj.cor_2pi = 1;
           end
        end
    end
    
    
    
end