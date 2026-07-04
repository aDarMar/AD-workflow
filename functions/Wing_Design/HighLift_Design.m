function [wing3Ddata] = HighLift_Design(aero_des)
%HIGHLIFT_DESIGN Summary of this function goes here
%   Detailed explanation goes here
fig_HL = figure('Name','High Lift Aerodynamic Data');
disp('%%%%%%%%%%%%%%%%%%%%%%%% WING DESIGN MODULE %%%%%%%%%%%%%%%%%%%%%%%%');
disp('%%%%%%%%%%%%%%%%%%%% HIGH LIFT DESIGN SUBMODULE %%%%%%%%%%%%%%%%%%%%');
m = 1; QT = 2; wing3Ddata = ProfileClass.empty
while QT > 1
    disp('Choose Flap and Slat Deflections')
    disp('1. Flap');
    deltaFs = input('>>');
    disp('2. Slat');
    deltaSs = input('>>');
    %% DATA Calculation
    disp(['Calculations will be performed at Mach = ',num2str( aero_des.low_speed.prf3DClean.M ) ]);
    if deltaSs == 0
        wing3Ddata(m)          = aero_des.low_speed.aero3Dwing( 'Take-Off',aero_des.low_speed.prf3DClean.M,deltaFs,deltaSs );
    else
        wing3Ddata(m)          = aero_des.low_speed.aero3Dwing( 'Landing',aero_des.low_speed.prf3DClean.M,deltaFs,deltaSs );
    end
    
    %% Results Display
    ax_HL(1) = subplot(1,2,1,'Parent',fig_HL); hold( ax_HL(1),'on' )
    ax_HL(2) = subplot(1,2,2,'Parent',fig_HL);
    a_min = min( aero_des.low_speed.prf3DClean.alpha0l,wing3Ddata(m).alpha0l );
    a_max = max( aero_des.low_speed.prf3DClean.alphamax,wing3Ddata(m).alphamax );
    alpha    = linspace( a_min,a_max,15 );
    CL_clean = aero_des.low_speed.lift_eval( alpha,aero_des.low_speed.prf3DClean );
    CL_hl    = aero_des.low_speed.lift_eval( alpha,wing3Ddata(m) );
    plot( ax_HL(1), alpha,CL_clean ); plot( ax_HL(1), alpha,CL_hl ); 
    disp('Do you want to save the results?')
    disp('1. Yes')
    disp('2. No');
    intr = input('>>');
    if intr == 1 
        m = m + 1;
    end
    disp('Quit?')
    disp('1. Yes')
    disp('2. No');
    intr = input('>>');
    if intr == 1
        QT = -1;
    else
        QT = 2;
    end
end


end

