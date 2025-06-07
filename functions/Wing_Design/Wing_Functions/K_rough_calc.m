function [k_rough] = K_rough_calc(Type_of_surface)
if Type_of_surface==1
    k_rough=0.0; % mm
elseif Type_of_surface==2
    k_rough=0.00127;% mm
elseif Type_of_surface==3
    k_rough=0.00406;% mm
elseif Type_of_surface==4
    k_rough=0.00635;% mm
elseif Type_of_surface==5
    k_rough=0.01016;% mm
end
end