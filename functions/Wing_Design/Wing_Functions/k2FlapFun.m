function [k2] = k2FlapFun( delta_f, flap_type, n_sections)

% cf_c DEVE essere un vettore con n_sections elementi
if length(delta_f) ~= n_sections
    error('delta_f deve essere un vettore con %d elementi.', n_sections);
end
% Flap type, inserisci: A per 1-slot
%                       B per 2-slots
%                       C per Fowler
%                       D per Split or plain


load ('Grafici\K2.mat');


switch flap_type
    case '1-slotted'
        k2=interp1(delta_flap_1_slot, K2_1_slot, delta_f);
        % Verifichiamo se gli input sono fuori dall'intervallo dei sample points
        % ed eventualmente li settiamo ai valori di bordo, in modo da avere come
        % valore di estrapolazione un valore pari a quello di bordo
        for i=1:n_sections
            if isnan(k2(i))
                if delta_f(i)>max(delta_flap_1_slot(:))
                    delta_f(i)=max(delta_flap_1_slot(:));
                end
                if delta_f(i)<min(delta_flap_1_slot(:))
                    delta_f(i)=min(delta_flap_1_slot(:));
                end
                warning('delta_f(%d) fuori intervallo: impostato al valore di bordo.', i);
                k2(i)=interp1(delta_flap_1_slot, K2_1_slot, delta_f(i));
            end
        end

    case '2-slotted'
        k2=interp1(delta_flap_2_slot, K2_2_slot, delta_f);
        for i=1:n_sections
            if isnan(k2(i))
                if delta_f(i)>max(delta_flap_2_slot(:))
                    delta_f(i)=max(delta_flap_2_slot(:));
                end
                if delta_f(i)<min(delta_flap_2_slot(:))
                    delta_f(i)=min(delta_flap_2_slot(:));
                end
                warning('delta_f(%d) fuori intervallo: impostato al valore di bordo.', i);
                k2(i)=interp1(delta_flap_2_slot, K2_2_slot, delta_f(i));
            end
        end
    case 'fowler'
        k2=interp1(delta_flap_fowler, K2_fowler, delta_f);
        for i=1:n_sections
            if isnan(k2(i))
                if delta_f(i)>max(delta_flap_fowler(:))
                    delta_f(i)=max(delta_flap_fowler(:));
                end
                if delta_f(i)<min(delta_flap_fowler(:))
                    delta_f(i)=min(delta_flap_fowler(:));
                end
                warning('delta_f(%d) fuori intervallo: impostato al valore di bordo.', i);
                k2(i)=interp1(delta_flap_fowler, K2_fowler, delta_f(i));
            end
        end
    case 'plain'
        k2=interp1(delta_flap_Split_Plain, K2_Split_Plain, delta_f);
        for i=1:n_sections
            if isnan(k2(i))
                if delta_f(i)>max(delta_flap_Split_Plain(:))
                    delta_f(i)=max(delta_flap_Split_Plain(:));
                end
                if delta_f(i)<min(delta_flap_Split_Plain(:))
                    delta_f(i)=min(delta_flap_Split_Plain(:));
                end
                warning('delta_f(%d) fuori intervallo: impostato al valore di bordo.', i);
                k2(i)=interp1(delta_flap_Split_Plain, K2_Split_Plain, delta_f(i));
            end
        end
    otherwise
        error('flap_type non valido. Usa: A, B, C oppure D');
end