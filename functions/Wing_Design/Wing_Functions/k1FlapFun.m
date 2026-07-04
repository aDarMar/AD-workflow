function [k1] = k1FlapFun( cf_c, flap_type, n_sections)

% cf_c DEVE essere un vettore con n_sections elementi
if length(cf_c) ~= n_sections
    error('cf_c deve essere un vettore con %d elementi.', n_sections);
end
% Flap type, inserisci: A per split, plain and 1-slot
%                       B per 2-slots or Fowler
load ('Grafici\K1.mat');

cf_c=100*cf_c;

if isequal('Flowler',flap_type) || isequal('2-slotted',flap_type)
    flap_type = 'fowler';
elseif isequal('1-slotted',flap_type) || isequal('plain',flap_type) || isequal('split',flap_type)
    flap_type = 'plain';
end

switch flap_type
    case 'plain'
        k1=interp1(cf_c_1_split_plain, K1__1_split_plain, cf_c); 
% Verifichiamo se gli input sono fuori dall'intervallo dei sample points
% ed eventualmente li settiamo ai valori di bordo, in modo da avere come
% valore di estrapolazione un valore pari a quello di bordo         
        for i=1:n_sections
            if isnan(k1(i))
               if cf_c(i)>max(cf_c_1_split_plain(:)) 
                   cf_c(i)=max(cf_c_1_split_plain(:));
                   %warning('cf_c(%d) fuori intervallo: impostato al valore di bordo.', i);
               end
               if cf_c(i)<min(cf_c_1_split_plain(:))
                   cf_c(i)=min(cf_c_1_split_plain(:)); 
                   %warning('cf_c(%d) fuori intervallo: impostato al valore di bordo.', i);
               end
               k1(i)=interp1(cf_c_1_split_plain, K1__1_split_plain, cf_c(i));
            end
        end
    case 'fowler'
        k1=interp1(cf_c_2_slot_fowler, K1_2_slot_fowler, cf_c); 
% Verifichiamo se gli input sono fuori dall'intervallo dei sample points
% ed eventualmente li settiamo ai valori di bordo, in modo da avere come
% valore di estrapolazione un valore pari a quello di bordo         
        for i=1:n_sections
            if isnan(k1(i))
               if cf_c(i)>max(cf_c_2_slot_fowler(:)) 
                   cf_c(i)=max(cf_c_2_slot_fowler(:));
                   %warning('cf_c(%d) fuori intervallo: impostato al valore di bordo.', i);
               end
               if cf_c(i)<min(cf_c_2_slot_fowler(:))
                   cf_c(i)=min(cf_c_2_slot_fowler(:)); 
                   %warning('cf_c(%d) fuori intervallo: impostato al valore di bordo.', i);
               end
               k1(i)=interp1(cf_c_2_slot_fowler, K1_2_slot_fowler, cf_c(i));
            end
        end
     otherwise
         error('flap_type non valido. Usa: A oppure B');
end