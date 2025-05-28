function [cf_turb] = cf_turb_calc(Re,Re_cut_off,M)
if Re>Re_cut_off
   cf_turb=0.455/( ((log10(Re_cut_off))^(2.58))*((1+0.144*M^2)^(0.65)));
else
   cf_turb=0.455/( ((log10(Re))^(2.58))*((1+0.144*M^2)^(0.65))); 
end