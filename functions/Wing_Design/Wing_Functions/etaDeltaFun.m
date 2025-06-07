function [etaDelta] = etaDeltaFun( delta_s, n_sections)

% cf_c DEVE essere un vettore con n_sections elementi
if length(delta_s) ~= n_sections
    error('delta_s deve essere un vettore con %d elementi.', n_sections);
end

load ('Grafici\ds_eta_delta.mat');

etaDelta=interp1(ds, eta_delta, delta_s);
% Verifichiamo se gli input sono fuori dall'intervallo dei sample points
% ed eventualmente li settiamo ai valori di bordo, in modo da avere come
% valore di estrapolazione un valore pari a quello di bordo
for i=1:n_sections
    if isnan(etaDelta(i))

        if delta_s(i)>max(ds(:))
            delta_s(i)=max(ds(:));
        end
        if delta_s(i)<min(ds(:))
            delta_s(i)=min(ds(:));
        end
        warning('delta_s(%d) fuori intervallo: impostato al valore di bordo.', i);
        etaDelta(i)=interp1(ds, eta_delta, delta_s(i));
    end
end
end