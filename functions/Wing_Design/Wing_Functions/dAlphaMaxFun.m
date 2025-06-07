function [dAlfa_max] = dAlphaMaxFun(Freccia_le, Dy_average)


load('grafici\m_Delta_alpha_max_chart_19.mat');
load("grafici\m_DeltaY_c_chart_5_19.mat");
load("grafici\m_Sweep_chart_5.19.mat");

% dAlfa_max1=interp2(m_Sweep,m_DeltaY_c,m_Delta_alpha_max,Freccia_le,Dy_average);
I = scatteredInterpolant(m_Sweep(:),m_DeltaY_c(:),m_Delta_alpha_max(:),...
    'linear', 'boundary');
dAlfa_max = I(Freccia_le,Dy_average);
%abs(1-dAlfa_max/dAlfa_max1)
end