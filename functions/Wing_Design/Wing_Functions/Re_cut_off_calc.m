function [Re_cut_off] = Re_cut_off_calc(l,k,M)
if (M<0.9)
Re_cut_off=38.21*(l/k)^(1.053);
elseif (M>=0.9)
Re_cut_off=44.62*((l/k)^(1.053))*M^(1.16);
end
end