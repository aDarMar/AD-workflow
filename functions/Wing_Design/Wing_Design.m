function [aero_des,outputArg2] = Wing_Design( aero_des )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

disp('%%%%%%%%%%%%%%%%%%%% WING DESIGN MODULE %%%%%%%%%%%%%%%%%%%%');
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
% Section Geometry
yob  = fscanf(f_id,'%f '); fgetl(f_id);
eps  = fscanf(f_id,'%f '); fgetl(f_id); % gets the first row
toc  = fscanf(f_id,'%f '); fgetl(f_id);
n_inp_geom = 9;
geom_vec = nan(3,n_inp_geom);
for i = 5:n_inp_geom
    %[2y/b,c,eps,t/c,] LER/c,x@max(t/c),xtr_Up,xtr_Low,dY
    geom_vec(:,i) = fscanf(f_id,'%f '); fgetl(f_id);
end
temp = fgetl(f_id);
% Secions Aerodynamics
tag = 'Aerodynamics';
if ~strcmp(temp,tag)
    error(['Error Reading ',tag])
end
ka   = fscanf(f_id,'%f '); fgetl(f_id);
temp = fgetl(f_id);
%Low Speed Data
tag = 'Low Speed';
if ~strcmp(temp,tag)
    error(['Error Reading ',tag, 'Aerodynamic'])
end
n_inp_aero = 11;
aero_vec_low = nan(3,n_inp_aero);
for i = 1:n_inp_aero
    % Mach,cla,cl0,cl*,clmax,alphamax,alpha0l,alpha*,cm_ac
    aero_vec_low(:,i) = fscanf(f_id,'%f '); fgetl(f_id);
end

% Low Speed Drag
tag = 'Low Speed Drag'; err_read(f_id,tag);
fgetl(f_id); low_speed_drag = iter_read(f_id,4);

%High Speed Data
tag = 'High Speed'; err_read(f_id,tag);
n_inp_aero = 11;
aero_vec_high = nan(3,n_inp_aero);
for i = 1:n_inp_aero
    % Mach,cla,cl0,cl*,clmax,alphamax,alpha0l,alpha*,cm_ac
    aero_vec_high(:,i) = fscanf(f_id,'%f '); fgetl(f_id);
end

% High Speed Drag
tag = 'High Speed Drag'; err_read(f_id,tag);
fgetl(f_id)
high_speed_drag = iter_read(f_id,4);

% High Lift Devices
tag = 'High-Lift Devices'; err_read(f_id,tag);
tag = 'Flaps'; err_read(f_id,tag); 
fgetl(f_id);  flaps_data = iter_read(f_id,4);

tag = 'Slats'; err_read(f_id,tag); %fgetl(f_id); 
fgetl(f_id); slats_data = iter_read(f_id,4);
fgetl(f_id); % Skips the row after

fclose(f_id);

%% Plantform Definition
[f_wing,ax_wing,yroot,ykink,ytip,croot,ckink,ctip,xLE_root,xLE_kink, ...
    xLE_tip,A1_Sw] = Wing_Plantform_fun( aero_des.Sw,aero_des.bw,aero_des.sweepw,aero_des.TRiw,...
    yob(2),aero_des.TRw,aero_des.AioSw );
aero_des.AioSw = A1_Sw; % Updating A1/Sw
%% Equivalent Wing Definition
aero_des = aero_des.equivalent_wing_def(ctip);
% hold on
% p1 = plot(ax_wing,[aero_des.equiv_wing.panels.root.xglob,aero_des.equiv_wing.panels.tip.xglob],...
%     [aero_des.equiv_wing.panels.root.yglob,aero_des.equiv_wing.panels.tip.yglob]); 
% p1 = plot(ax_wing,[aero_des.equiv_wing.panels.root.xglob,aero_des.equiv_wing.panels.tip.xglob]+...
%     [aero_des.equiv_wing.panels.root.c,aero_des.equiv_wing.panels.tip.c],...
%     [aero_des.equiv_wing.panels.root.yglob,aero_des.equiv_wing.panels.tip.yglob]);
% p1 = plot(ax_wing,[aero_des.equiv_wing.panels.tip.xglob,aero_des.equiv_wing.panels.tip.xglob]...
%     +[0,aero_des.equiv_wing.panels.tip.c],...
%     [aero_des.equiv_wing.panels.root.yglob,aero_des.equiv_wing.panels.tip.yglob]);
% hold off
% Building the Geometric Vector
%    2y/b,c,eps,t/c,LER/c,x@max(t/c),xtr_Up,xtr_Low,dY
i = 1;
geom_vec(1,i) = yroot; geom_vec(2,1) = ykink; geom_vec(3,i) = ytip;
i = 2;
geom_vec(1,i) = croot; geom_vec(2,i)   = ckink; geom_vec(3,i) = ctip;
i = 3; geom_vec(:,i) = eps(:)'; 

%% Airfoil Selection
toc = Airfoil_Selection( aero_des,0.02,ka(1),toc,...
    geom_vec(:,2),geom_vec(:,1) ); % Temporary ka solution
i = 1; geom_vec(:,i) = geom_vec(:,i)./geom_vec(end,i);
i = 4; geom_vec(:,i) = toc(:)';


%% Wing Circulation
apexC = [ aero_des.wingapex.x,aero_des.wingapex.y,aero_des.wingapex.z ];
alpha_v = -3:18;
% Low Speed Wing Object Definition
m = 7; M = 7; % Spanwise section for Weissinger Method
HLflag = HL_build( flaps_data,slats_data );
aero_des.low_speed = PaneledWing( m,M,geom_vec,aero_vec_low,aero_des.bw,...
    aero_des.sweepw,aero_des.dihedralw,aero_des.iw,apexC,aero_vec_low(1,1),...
    HLflag,[flaps_data(:,2:4);slats_data(:,2:4)],nan,1 ); % first column of flap and slat data is the flap ID
% 3D data calculation and estimation
aero_des.low_speed.prf3DClean  = aero_des.low_speed.aero3Dwing( 'clean', aero_des.low_speed.panels(1).root.M );
% High Lift Definition
aero_des.low_speed.wing3Ddata  = HighLift_Design(aero_des);
% aero_des.low_speed.wing3Ddata(m)          = aero_des.low_speed.aero3Dwing( 'Take-off',M(m),deltaFs(m),deltaSs(m) );
% aero_des.low_speed.wing3Ddata(m)          = aero_des.low_speed.aero3Dwing( 'Landing',M(m),deltaFs(m),deltaSs(m) );

% Drag Definition
aero_des.low_speed.meanprofile.poly_drag  = aero_des.low_speed.poly_drag( low_speed_drag(:,1),low_speed_drag(:,2:4) ); % Defines the interpolating function for cd avg
CDl = aero_des.CDlow_Mach( alpha_v );
%cm  = aero_des.low_speed.cm_alpha( alpha_v );
% High Speed
aero_des.high_speed = PaneledWing( m,M,geom_vec,aero_vec_high,aero_des.bw,...
    aero_des.sweepw,aero_des.dihedralw,aero_des.iw,apexC,aero_vec_high(1,1) );
% 3D data calculation and estimation
aero_des.high_speed.prf3DClean            = aero_des.high_speed.aero3Dwing( 'clean', aero_des.high_speed.panels(1).root.M );
aero_des.high_speed.meanprofile.poly_drag = aero_des.high_speed.poly_drag( high_speed_drag(:,1),high_speed_drag(:,2:4) );
CDh = aero_des.CDtransonic( alpha_v,aero_des.TLARs.cruise.M );
aero_des.Mdd;

CL_low  = aero_des.low_speed.lift_eval(alpha_v,aero_des.low_speed.prf3DClean);
CL_high = aero_des.high_speed.lift_eval(alpha_v,aero_des.high_speed.prf3DClean );
aero_des.plot_fun(alpha_v,CL_low,CL_high,CDl, CDh)
aero_des.low_speed.wing_circ(3)
end

function out = iter_read(f_id,n_i)
    %iter_read: function that reads rows of data from files.
    %INPUT:
    %   f_id: ID of file we want to reaad from
    %   n_i: size of data we expect to read from a column
    temp = zeros(1,n_i); i = 1;
    temp = fscanf(f_id,'%f'); fgetl(f_id);
    while ~isempty(temp)
        out(i,:) = temp(:)'; 
        temp = fscanf(f_id,'%f '); fgetl(f_id);
        i = i+1;
    end
    %out   = out(2:end,:); % Excludes firt row of all zeros
end

function err_read(f_id,tag)
    temp = fgetl(f_id);
    if ~strcmp(temp,tag)
        error(['Error Reading ', tag ,' Data'])
    end
end

function HLF = HL_build( flaps_data,slats_data )
    n_fl = length( flaps_data(:,1) ); 
    Qo = floor(n_fl / 2);  % Quoziente intero
    R = mod(n_fl, 2);    % Resto
    if R > 0
        warning('Number of Flap sections is odd');
    end
    HLF = repmat( {'flaps'},Qo,1);
    n_sl = length( slats_data(:,1) ); 
    Q = floor(n_sl / 2);  % Quoziente intero
    R = mod(n_sl, 2);    % Resto
    if R > 0
        warning('Number of Slat sections is odd');
    end
    for i = 1:Q
     HLF{Qo+i} = 'slats';   
    end
    
    
end