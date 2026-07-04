function [Kc] = kcFun(avg2Dalfad, ARw, n_flap)

% twoYob deve essere un vettore con 2 elementi
if length(avg2Dalfad) ~= n_flap
    error('avg2Dalfad deve essere un vettore con 2 elementi, ordinati inboard-outboard.');
end

load ('Grafici\m_AR_K_C_chart.mat');
load ('Grafici\m_Alpha_delta_k_c_CHART.mat'); % ERRORE NEL FILE DI INPUT, HA UNA COLONNA DUPLICATA
load ('Grafici\m_KC_K_C_chart.mat');

I = scatteredInterpolant(m_AR(:), m_Alpha_delta(:), m_KC(:),'linear','boundary');
Kc = I(avg2Dalfad, ARw);

%kc= interp2(m_AR, m_Alpha_delta, m_KC, avg2Dalfad, ARw);

end