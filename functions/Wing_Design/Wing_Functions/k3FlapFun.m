function [k3] = k3FlapFun( df_dfref, flap_type, n_sections)
%k3FlspFun calcola il coefficiente K3

% cf_c DEVE essere un vettore con n_sections elementi
if length(df_dfref) ~= n_sections
    error('df_dfref deve essere un vettore con %d elementi.', n_sections);
end
% Flap type, inserisci: A per 1-slot or Fowler
%                       B per 2-slots


load ('Grafici\K3.mat');

if isequal('fowler',flap_type) || isequal('1-slotted',flap_type)
    flap_type = 'fowler';
end


switch flap_type

    case 'fowler'
        k3=interp1(df_df_ref_1_slot, K3_1_slot, df_dfref);
        % Verifichiamo se gli input sono fuori dall'intervallo dei sample points
        % ed eventualmente li settiamo ai valori di bordo, in modo da avere come
        % valore di estrapolazione un valore pari a quello di bordo
        for i=1:n_sections
            if isnan(k3(i))
                if df_dfref(i)>max(df_df_ref_1_slot(:))
                    df_dfref(i)=max(df_df_ref_1_slot(:));
                end
                if df_dfref(i)<min(df_df_ref_1_slot(:))
                    df_dfref(i)=min(df_df_ref_1_slot(:));
                end
                warning('df_dfref(%d) fuori intervallo: impostato al valore di bordo.', i);
                k3(i)=interp1(df_df_ref_1_slot, K3_1_slot, df_dfref(i));
            end
        end

    case '2-slotted'
        k3=interp1(df_df_ref_2_slot, K3__2_slot, df_dfref);
        for i=1:n_sections
            if isnan(k3(i))
                if df_dfref(i)>max(df_df_ref_2_slot(:))
                    df_dfref(i)=max(df_df_ref_2_slot(:));
                end
                if df_dfref(i)<min(df_df_ref_2_slot(:))
                    df_dfref(i)=min(df_df_ref_2_slot(:));
                end
                warning('df_dfref(%d) fuori intervallo: impostato al valore di bordo.', i);
                k3(i)=interp1(df_df_ref_2_slot, K3__2_slot, df_dfref(i));
            end
        end
    otherwise
        error('flap_type non valido. Usa: A oppure B');
end