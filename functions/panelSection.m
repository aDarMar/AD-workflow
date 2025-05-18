classdef panelSection < ProfileClass
   properties
       eta
       phi
   end
    methods
        function obj = panelSection(geomV,aeroV,M,phi,b)
           obj@ProfileClass(geomV,aeroV,M) % call superclass constructor
           obj.phi   = phi;
           obj.eta   = cos(obj.phi);
           obj.yglob = obj.eta*b*0.5;
        end
    end
    
    
    
end