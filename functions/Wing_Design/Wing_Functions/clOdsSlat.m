function [dCloDs] = clOdsSlat(csoc)
%clOdsSlat Valori della derivata Cl/ds rispetto il rapporto delle corde tra
%slat e profilo. 
%   Rate of change of airfoil lift coefficient with slat 
% chord ratio
% csoc: rapporto corda sezione con slat e senza

load('Grafici\slat_cld_ds.mat');

if csoc > max(cs_c(:))
    csoc = max(cs_c(:));
end
if csoc < min(cs_c(:))
    csoc = min(cs_c(:));
end

dCloDs = interp1(cs_c,cld_ds,csoc,'linear');


end