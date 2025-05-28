function [etaMax] = etaMaxFun( LER_t, n_sections)

% cf_c DEVE essere un vettore con n_sections elementi
if length(LER_t) ~= n_sections
    error('LER_t deve essere un vettore con %d elementi.', n_sections);
end
load ('Grafici\eta_max_LER.mat');
etaMax=interp1(LER__t_c, eta_max, LER_t);
% Verifichiamo se gli input sono fuori dall'intervallo dei sample points
% ed eventualmente li settiamo ai valori di bordo, in modo da avere come
% valore di estrapolazione un valore pari a quello di bordo
for i=1:n_sections
    if isnan(etaMax(i))
        if LER_t(i)>max(LER__t_c(:))
            LER_t(i)=max(LER__t_c(:));
        end
        if LER_t(i)<min(LER__t_c(:))
            LER_t(i)=min(LER__t_c(:));
        end
        %warning('LER_t(%d) fuori intervallo: impostato al valore di bordo.', i);
        etaMax(i)=interp1(LER__t_c, eta_max, LER_t(i));
    end
end

end