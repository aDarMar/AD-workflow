function [CL_max] = CL_max_fun(Mean_Cl_max,sweep,dy,~, M_inf)
%CL_max_fun: ricava il CL max per un'ala tridimensionale in configurazione
%pulita.
    dy              = dy*100;
    Delta_CL_max    = dCLMaxFun( sweep, dy, M_inf );
    CL_max_o_Cl_max = clCLFun(sweep,dy);
    CL_max          = ( CL_max_o_Cl_max*Mean_Cl_max ) + Delta_CL_max;

end
