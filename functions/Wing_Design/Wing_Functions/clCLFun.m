function [CLocl] = clCLFun(sweep,dy,c)
load ('grafici\m_CL_max_cl_max_chart_5_17.mat');
load ('grafici\m_DeltaY_chart_5_17.mat');
load ('grafici\m_Sweep_chart_5_17.mat');

dy_c=dy/c;
I = scatteredInterpolant(m_Sweep(:),m_DeltaY(:),m_CL_max_Cl_max(:),'linear','boundary');
CLocl = I(sweep,dy_c);

end