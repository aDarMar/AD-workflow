close all; clear; clc

addpath('functions')
addpath('functions\Wing_Design')
addpath('functions\Wing_Design\Wing_Functions')
addpath('functions\Wing_Design\Wing_Functions\grafici')
addpath('functions\Weight_estimation');
addpath('functions\Sizing_Plot');
addpath('functions\Classes');
addpath('functions\Wing_Design\u_TR');
M = 7; m_in = 2; m = 7;
geom_vec = nan(m_in,9); aero_vec = nan(m_in,9);
aero_vec(:,1) = 0; aero_vec(:,2) = 2*pi;
%temp = [ 3.152240935; 3.585786438; 4.234633135;5 ];
%geom_vec(:,1) = [ temp;flip(temp(1:end-1)) ];
geom_vec(:,1) = [1,0]; geom_vec(:,2) = [3,5];


% OLD
%wing_d = aero_influence_coeffs( m,M,geom_vec,aero_vec,b,sweep );
%wing_d.basicLoad_coeffs
% NEW
CHS = 5;
switch CHS
    case 1
        b = 36; sweep = 20; dihedral = 4; iang = 0; apexC = [0,0,0]; Mach = 0.1;
        n_in = 5;
        geom_vec = nan( n_in,9 ); aero_vec = nan( n_in,9 );
        geom_vec( :,1 ) = [0,0.111111111,0.444444444,0.777777778,1]';
        geom_vec( :,2 ) = [6,5.272059531,3.969873029,2.844013508,2.093440494]';
    case 2
        n_in = 6;
        geom_vec = nan( n_in,9 ); aero_vec = nan( n_in,9 );
        geom_vec( :,1 ) = [0,2,8,10,14,16]';
        geom_vec( :,2 ) = [5,4.75,4,3.75,3.25,3]';
        b = geom_vec( end,1 )*2; sweep = 21.5649; dihedral = 4; iang = 0; apexC = [0,0,0]; Mach = 0.1;
        geom_vec( :,1 ) = geom_vec( :,1 )/geom_vec( end,1 );
    case 3
        n_in = 5;
        geom_vec = nan( n_in,9 ); aero_vec = zeros( n_in,9 );
        geom_vec( :,1 ) = [0,2,3,6,9]';
        geom_vec( :,2 ) = 3.35/(2*pi)*[4,3.555555556,3.333333333,2.666666667,2]';
        geom_vec( :,3 ) = [0,-0.222222222,-0.333333333,-0.666666667,-1]';
        b = geom_vec( end,1 )*2; sweep = 46.5482; dihedral = 4; iang = 0; apexC = [0,0,0]; Mach = 0.1;
        geom_vec( :,1 ) = geom_vec( :,1 )/geom_vec( end,1 );
    case 4
        n_in = 5;
        geom_vec = nan( n_in,9 ); aero_vec = zeros( n_in,9 );
        geom_vec( :,1 ) = [0,2,3,6,9]';
        geom_vec( :,2 ) = [4,3.555555556,3.333333333,2.666666667,2]';
        geom_vec( :,3 ) = [0,-0.222222222,-0.333333333,-0.666666667,-1]';
        aero_vec( :,2 ) = (3.35/6)*[6,6,6,6,6]';%[6,5.9,5.9,5.9,5.9]';
        b = geom_vec( end,1 )*2; sweep = 46.5482; dihedral = 4; iang = 0; apexC = [0,0,0]; Mach = 0.1;
        geom_vec( :,1 ) = geom_vec( :,1 )/geom_vec( end,1 );
    case 5
        n_in = 5;
        geom_vec = nan( n_in,9 ); aero_vec = zeros( n_in,9 );
        geom_vec( :,1 ) = [0,2,3,6,9]'; 
        geom_vec( :,2 ) = [4,3.555555556,3.333333333,2.666666667,2]';
        geom_vec( :,3 ) = [0,-0.222222222,-0.333333333,-0.666666667,-1]';
        aero_vec( :,1 ) = (sqrt(0.19)/6)*[6,6,6,6,6]';
        aero_vec( :,2 ) = (5.28/6)*[6,6,6,6,6]';
        b = geom_vec( end,1 )*2; sweep = 43.698; dihedral = 4; iang = 0; apexC = [0,0,0]; Mach = aero_vec( 1,2 );
        geom_vec( :,1 ) = geom_vec( :,1 )/geom_vec( end,1 );
        
end
des_wing = PaneledWing( m,M,geom_vec,aero_vec,b,sweep,dihedral,iang,apexC,Mach );
%des_wing = des_wing.aeroDef;

%% Wing_Design Debug
main_fold = cd;
TLARs_path = [main_fold,'\tlars\TLARs.txt'];
Des_Air = Air_Design( TLARs_path,main_fold );

Des_Air.bw = 34.67; Des_Air.Sw = 133.54;
Des_Air.TLARs.cruise.M = 0.785;
Des_Air.CL_cr = 0.32; Des_Air.CLmax_cr = 1.50;
Wing_Design( Des_Air );
