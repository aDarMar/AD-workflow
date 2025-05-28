function [etadelta_f] = etadelta_fFun(delta_f, flap_type, n_sections)

% Gli input devono essere vettori con n_sections elementi
if length(delta_f) ~= n_sections
    error('delta_f deve essere un vettore con %d elementi.', n_sections);
end

load ('Grafici\eta_delta_f.mat');

switch flap_type
    case '1-slotted'
        if delta_f > max(delta_flap_1_slot)
            delta_f = max(delta_flap_1_slot);
        end
        if delta_f < min(delta_flap_1_slot)
            delta_f = min(delta_flap_1_slot);
        end
        etadelta_f = interp1(delta_flap_1_slot, eta_delta_1_slot,delta_f,'linear');
    case '2-slotted'
        if delta_f > max(delta_flap_2_slot)
            delta_f = max(delta_flap_2_slot);
        end
        if delta_f < min(delta_flap_2_slot)
            delta_f = min(delta_flap_2_slot);
        end
        etadelta_f = interp1(delta_flap_2_slot, eta_delta_2_slot,delta_f,'linear');
    case '3-slotted'
        if delta_f > max(delta_flap_3_slot)
            delta_f = max(delta_flap_3_slot);
        end
        if delta_f < min(delta_flap_3_slot)
            delta_f = min(delta_flap_3_slot);
        end
        etadelta_f = interp1(delta_flap_3_slot, eta_delta_3_slot,delta_f,'linear');
    case 'fowler'
        if delta_f > max(delta_flap_fowler)
            delta_f = max(delta_flap_fowler);
        end
        if delta_f < min(delta_flap_fowler)
            delta_f = min(delta_flap_fowler);
        end
        etadelta_f = interp1(delta_flap_fowler,eta_delta_fowler,delta_f,'linear');
    otherwise
        error("Flap non Riconosciuto")
end

%etadelta_f=interp1(delta_flap_1_slot, eta_delta_1_slot, delta_f);