function [etad] = etaDeltaFlapFun(df,flaptype)
%etaDeltaFlapFun interpola per ottenere i punti del grafico slide 16
load('Grafici\Deltac_cf_delta_f.mat')

if isequal('3-slotted',flaptype)
    flaptype = '2-slotted';
end

switch flaptype
    case 'fowler'
        if df> max(delta_flap_fowler)
            df = max(delta_flap_fowler);
        end
        if df<min(delta_flap_fowler)
            df = min(delta_flap_fowler);
        end
        etad = interp1(delta_flap_fowler,Delta_cf_c_fowler,df,"linear");
    case '1-slotted'
        if df> max(delta_flap_1_slot)
            df = max(delta_flap_1_slot);
        end
        if df<min(delta_flap_1_slot)
            df = min(delta_flap_1_slot);
        end
        etad = interp1(delta_flap_1_slot,Delta_cf_c_1_slot,df,"linear");

    case '2-slotted'
        if df> max(delta_flap_2_slot)
            df = max(delta_flap_2_slot);
        end
        if df<min(delta_flap_2_slot)
            df = min(delta_flap_2_slot);
        end
        etad = interp1(delta_flap_2_slot,Delta_cf_c_2_slot,df,"linear");
    case 'plain'
        if df> max(delta_flap_plain_slot)
            df = max(delta_flap_plain_slot);
        end
        if df<min(delta_flap_plain_slot)
            df = min(delta_flap_plain_slot);
        end
        etad = interp1(delta_flap_plain_slot,Delta_cf_cplain_slot,df,"linear");
    otherwise
        error("Tipo di flap non riconosciuto")
end
end