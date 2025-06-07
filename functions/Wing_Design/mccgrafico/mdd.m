function MDDad(tc_mean, sweep_LE, CL_cruise,sweep_c4 )
% Trova il MDD con l'approccio che usa ADAS, ovvero Stanford-Kroo
%Input:
%tc_mean = t/c del profilo medio;
%sweep_LE = angolo di freccia al Leading Edge
%CL_cruise = CL del profilo 3D in crociera
%sweep_c4 = angolo di freccia a c/4
fatt = CL_cruise/(cos(sweep_LE/57.3))^2;
tcperp = tc_mean/cos(sweep_LE/57.3);
%per caricare il grafico digitalizzato
files =dir('*.csv');
X = [];
Y = [];
Z = [];
for k = 1:length(files)
    cl = files(k).name;
    data = load(cl);
    X = [X;data(:,1)];
    Y= [Y;data(:,2)];
    Z =  [Z;data(:,3)];
end 

F = scatteredInterpolant(X,Z,Y);
y = F(tcperp,fatt);
MCC = (y/cos(sweep_LE))+0.06;

MDDad = MCC*[1.02+0.08*(1-cos(1-cos(sweep_c4)))];
end
