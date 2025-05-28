function [CL_max] = CL_max_fun(Mean_Cl_max,sweep,dy,c, M_inf)
%CL_max_fun: ricava il CL max per un'ala tridimensionale in configurazione
%pulita.
    Delta_CL_max = dCLMaxFun(sweep, dy,M_inf);
    CL_max_o_Cl_max = clCLFun(sweep,dy,c);
    CL_max=(CL_max_o_Cl_max*Mean_Cl_max)+Delta_CL_max;

end
