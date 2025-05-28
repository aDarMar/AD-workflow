function [dClmaxbase] = dClMaxBaseFun(t_c, flap_type, n_sections)

% t_c DEVE essere un vettore con n_sections elementi
if length(t_c) ~= n_sections
    error('t_c deve essere un vettore con %d elementi.', n_sections);
end
% Flap type, inserisci: A per best 2-slots
%                       B per average 2-slots or Fowler
%                       C per NACA 2-slots or 1-slot
%                       D per split or plain

load ('Grafici\t_Delta_CL_max.mat');

t_c=100*t_c;

switch flap_type

    case '2-slotted'
        dClmaxbase=interp1(t_c_A, DeltaCLmaxbase_A, t_c);
        
% Verifichiamo se gli input sono fuori dall'intervallo dei sample points
% ed eventualmente li settiamo ai valori di bordo, in modo da avere come
% valore di estrapolazione un valore pari a quello di bordo
        for i=1:n_sections
               
             if isnan(dClmaxbase(i))

               if t_c(i)>max(t_c_A(:))
                   t_c(i)=max(t_c_A(:));
                   warning('t_c(%d) fuori intervallo: impostato al valore di bordo.', i);
               end
               if t_c(i)<min(t_c_A(:))
                   t_c(i)=min(t_c_A(:)); 
                   warning('t_c(%d) fuori intervallo: impostato al valore di bordo.', i);
               end
       
               dClmaxbase(i)=interp1(t_c_A, DeltaCLmaxbase_A, t_c(i));
            end
        end



    case 'fowler'
        dClmaxbase=interp1(t_c_B, DeltaCLmaxbase_B, t_c);
        for i=1:n_sections
               
 
            if isnan(dClmaxbase(i))

               if t_c(i)>max(t_c_B(:))
                   t_c(i)=max(t_c_B(:));
                    warning('t_c(%d) fuori intervallo: impostato al valore di bordo.', i);
               end
               if t_c(i)<min(t_c_B(:))
                   t_c(i)=min(t_c_B(:));
                    warning('t_c(%d) fuori intervallo: impostato al valore di bordo.', i);
               end
       
               dClmaxbase(i)=interp1(t_c_B, DeltaCLmaxbase_B, t_c(i));
            end
        end
        
    case 'C'
        dClmaxbase=interp1(t_c_C, DeltaCLmaxbase_C, t_c);
        for i=1:n_sections
               
             
            if isnan(dClmaxbase(i))

               if t_c(i)>max(t_c_C(:))
                   t_c(i)=max(t_c_C(:));
                   warning('t_c(%d) fuori intervallo: impostato al valore di bordo.', i);
               end
               if t_c(i)<min(t_c_C(:))
                   t_c(i)=min(t_c_C(:));
                   warning('t_c(%d) fuori intervallo: impostato al valore di bordo.', i);
               end
       
               dClmaxbase(i)=interp1(t_c_C, DeltaCLmaxbase_C, t_c(i));
            end
        end
        
    case 'D'
        dClmaxbase=interp1(t_c_D, DeltaCLmaxbase_D, t_c);
        for i=1:n_sections
               
 
            if isnan(dClmaxbase(i))

               if t_c(i)>max(t_c_D(:))
                   t_c(i)=max(t_c_D(:));
                   warning('t_c(%d) fuori intervallo: impostato al valore di bordo.', i);
               end
               if t_c(i)<min(t_c_D(:))
                   t_c(i)=min(t_c_D(:));
                   warning('t_c(%d) fuori intervallo: impostato al valore di bordo.', i);
               end
       
               dClmaxbase(i)=interp1(t_c_D, DeltaCLmaxbase_D, t_c(i));
            end
        end
    otherwise
         error('flap_type non valido. Usa: A, B, C o D.');
end