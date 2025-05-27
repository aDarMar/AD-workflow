function [aero_des,outputArg2] = Wing_Design( aero_des )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%% Reading Input Data
pth = [aero_des.main_fold,'\tlars\Geometry.txt'];
f_id = fopen(pth);
% First Line Check
temp = fgetl(f_id);
tag = 'Design Parameters';
if ~strcmp(temp,tag)
    error('File format inavalid')
end
tag = 'Sections';
while ~strcmp(temp,tag)
    temp = fgetl(f_id);
end
tag = 'Geometry';
while ~strcmp(temp,tag)
    temp = fgetl(f_id);
end
fgetl(f_id);
eps  = fscanf(f_id,'%f '); fgetl(f_id); % gets the first row 
toc  = fscanf(f_id,'%f '); fgetl(f_id);
yob  = fscanf(f_id,'%f '); fgetl(f_id);
fclose(f_id);

%% Plantform Definition
[yroot,ykink,ytip,croot,ckink,ctip,xLE_root,xLE_kink, ...
    xLE_tip] = Wing_Plantform_fun( aero_des.Sw,aero_des.bw,aero_des.sweepw,aero_des.TRiw,...
    yob(2),aero_des.TRw,aero_des.AioSw );
%% Equivalent Wing Definition
aero_des = aero_des.equivalent_wing_def(ctip);
%    2y/b,c,eps,t/c,LER/c,x@max(t/c),xtr_Up,xtr_Low,dY
% Building the Geometric Vector 
i = 1;
geom_vec(i,1) = yroot; geom_vec(i,2) = ykink; geom_vec(i,3) = ytip;
i = 2;
geom_vec(i,1) = croot; geom_vec(i,2)   = ckink; geom_vec(i,3) = ctip;
i = 3; geom_vec(i,:) = eps(:)'; 

%% Airfoil Selection
toc = Airfoil_Selection( aero_des.TLARs.cruise.M,aero_des.sweepw,...
    aero_des.CL_cr,aero_des.CLmax_cr,0.05,ka,toc,...
    geom_vec(:,2),geom_vec(:,1) );
i = 1; geom_vec(i,:) = geom_vec(i,:)./geom_vec(i,end);
i = 4; geom_vec(i,:) = toc(:)';

%% Wing Circulation
% Low Speed
m = 7; M = 7;
aero_des.low_speed = PaneledWing( m,M,geom_vec,aero_vec_low,aero_des.bw,...
    sweep,dihedral,iang,apexC,Mach );
% High Speed
aero_des.high_speed = PaneledWing( m,M,geom_vec,aero_vec_high,aero_des.bw,...
    sweep,dihedral,iang,apexC,Mach );
end