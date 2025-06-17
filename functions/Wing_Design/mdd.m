function [MDDad,MCC]= mdd(tc_mean, sweep_LE, CL_cruise,sweep_c4,main_fold)
% Trova il MDD con l'approccio che usa ADAS, ovvero Stanford-Kroo
%Input:
%tc_mean = t/c del profilo medio;
%sweep_LE = angolo di freccia al Leading Edge
%CL_cruise = CL del profilo 3D in crociera
%sweep_c4 = angolo di freccia a c/4
fatt = CL_cruise./(cos(sweep_LE/57.3)^2);fatt =fatt(:);
tcperp =ones(length(CL_cruise),1)*tc_mean/cos(sweep_LE/57.3);
files = [main_fold,'\functions\Wing_Design\mccgrafico\cl0..csv'];
            MAT = readmatrix(files);
            F = scatteredInterpolant( MAT(:,1),MAT(:,3),MAT(:,2) );
          F.ExtrapolationMethod = 'boundary';
y = F(tcperp,fatt);
MCC = (y./cos(sweep_LE/57.3))+0.06;

MDDad = MCC.*(1.02+0.08*(1-cos(1-cos(sweep_c4/57.3))));

end
