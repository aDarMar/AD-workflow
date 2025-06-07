function [Kb] = kbFun(twoYob, taper)

% twoYob deve essere un vettore con 2 elementi
if length(twoYob) ~= 2
    error('twoYob deve essere un vettore con 2 elementi, ordinati inboard-outboard.');
end

load ('Grafici\m_eta_flap_KB_chart.mat');
load ('Grafici\m_taper_ratio_KB_chart.mat');
load ('Grafici\m_KB_KB_chart.mat');

I = scatteredInterpolant(m_eta_flap(:), m_taper_ratio(:), m_KB(:),'linear','boundary');
kb = I(twoYob,ones(1,2)*taper);
% kb= interp2(m_eta_flap, m_taper_ratio, m_KB, twoYob, taper);

% Verifichiamo se gli input sono fuori dall'intervallo dei sample points
% ed eventualmente li settiamo ai valori di bordo, in modo da avere come
% valore di estrapolazione un valore pari a quello di bordo         
        % for i=1:2
        % 
        %     if isnan(kb(i))
        % 
        %        if twoYob(i)>max(m_eta_flap(:)) 
        %           twoYob(i)=max(m_eta_flap(:));
        %           warning('twoYob(%d) fuori intervallo: impostato al valore di bordo.', i); 
        %        end
        %        if twoYob(i)<min(m_eta_flap(:))
        %           twoYob(i)=min(m_eta_flap(:)); 
        %           warning('twoYob(%d) fuori intervallo: impostato al valore di bordo.', i); 
        %        end
        %        if taper<min(m_taper_ratio(:))
        %           taper=min(m_taper_ratio(:)); 
        %           warning('taper fuori intervallo: impostato al valore di bordo.', i); 
        %        end  
        %        if taper>max(m_taper_ratio(:))
        %           taper=max(m_taper_ratio(:)); 
        %           warning('taper fuori intervallo: impostato al valore di bordo.', i); 
        %        end  
        % 
        % 
        %        kb(i)= interp2(m_eta_flap, m_taper_ratio, m_KB, twoYob(i), taper);
        %     end
        % end

Kb=kb(2)-kb(1);

end