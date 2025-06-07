function  [ax1,ax2] = sizing_plot_Climb( ax1,ax2,...
    CLmaxes_TO,CLmaxes_LND,CLmaxes_CR,da_c,...
    CD0,TisaoT50 ,Roc_req)
%sizing_plot_Climb Function for plotting climb requirements
%INPUT
%   ax1: axis object for plotting limitations in [kg/m^2]
%   ax2: axis object for plotting limitations in [lb/ft^2]
%   CLmaxes_TO: vector of CL maxes for T/O
%   CLmaxes_LND: vector of CL maxes for T/O
%   CLmaxes_CR: vector of CL maxes for cruise
%   e: assumed Oswald factor
%   de_TO:  e decrement in T/O
%   de_LND: e decrement in LND
%   CD0: assumed drag coefficient for clean configuration
%   dCD0_f_TO: increment of drag coefficient due to flaps in T/O setting
%   dCD0_f_LND: increment of drag coefficient due to flaps in LND setting
%   dCD0_f_APP: increment of drag coefficient due to flaps in approach setting
%   dCD0_lgs: increment of drag coefficient due to extended landing gears
%   AR: assumed aspect ratio
%   neng: number of engines
%   WlandoMTOM: W_LND / W_TO
%   T0oTmo: max static thrust over max continous thrust ratio
%   TisaoT50: thrust reduction ratio between ISA condition and ISA+50
%   WcroWTo: weight fraction of W@end of cruise / W TO
%   WcroWTo_in: weight fraction of W@initial cruise / W TO
%   Roc_req: vector of max RoCs required at given height: [ RoC(i),h(i) ]
nCLs = length(CLmaxes_TO); 
%% Graphics
colors = [
    0.22, 0.49, 0.72;  % blu brillante
    0.53, 0.81, 0.92;  % azzurro chiaro
    0.27, 0.51, 0.71;  % blu acciaio
    0.00, 0.74, 0.83;  % azzurro-verde (cyan)
    0.00, 0.00, 0.55;  % blu scuro
    0.42, 0.35, 0.80;  % blu-grigio
    ];

if nCLs > length( colors )
    aux_c = random( nCLmaxes-colors,3 );
    colors = [colors;aux_c];
end

%% Conditions from Regulations:
% 1 - first segment of climb
CGD_norm(1)   = 0.012;              % Req. Climb Gradient
VoVmin(1)     = 1.2;                % Req V in terms of V/V_stall_TO
dCD0_flaps(1) = da_c.TLARs.dCD0_f_TO;          % CD0 increment due to flaps in T/O config
dCD0_LG(1)    = 0;                  % Land. Gear Retracted
e_vet(1)      = da_c.TLARs.e+da_c.TLARs.de_TO;
eng_w(1)      = (da_c.TLARs.nengine-1) / da_c.TLARs.nengine;    % OEI
WoWTO(1)      = 1;                  % MTOM
CL = CLmaxes_TO(:);                 % T/O 
% 2 - transitio to climb
CGD_norm(2)   = 0;                  % Req. Climb Gradient
VoVmin(2)     = 1.1;                %Req V in terms of V/V_stall_TO
dCD0_flaps(2) = da_c.TLARs.dCD0_f_TO;          % CD0 increment due to flaps in T/O config
dCD0_LG(2)    = da_c.TLARs.dCD0_lgs;           % Land. Gear Extended
e_vet(2)      = da_c.TLARs.e+da_c.TLARs.de_TO;
eng_w(2)      = (da_c.TLARs.nengine-1) / da_c.TLARs.nengine;    % OEI
WoWTO(2)      = 1;                  % MTOM
CL = [CL,CLmaxes_TO(:)];            % T/O 
% 3 - second segment
CGD_norm(3)   = 0.024;              % Req. Climb Gradient
VoVmin(3)     = 1.2;                % Req V in terms of V/V_stall_L
dCD0_flaps(3) = da_c.TLARs.dCD0_f_TO;          % CD0 increment due to flaps
dCD0_LG(3)    = 0;                  % Land. Gear Extended
e_vet(3)      = da_c.TLARs.e+da_c.TLARs.de_TO;
eng_w(3)      = (da_c.TLARs.nengine-1) / da_c.TLARs.nengine;    % OEI
WoWTO(3)      = 1;                  % MTOM
CL = [CL,CLmaxes_TO(:)]; % T/O 
% 4 - en-route Climb
CGD_norm(4)   = 0.012;              % Req. Climb Gradient
VoVmin(4)     = 1.25;               % Req V in terms of V/V_stall_L
dCD0_flaps(4) = 0;                  % CD0 increment due to flaps
dCD0_LG(4)    = 0;                  % Land. Gear Retracted
e_vet(4)      = da_c.TLARs.e;
eng_w(4)      = (da_c.TLARs.nengine-1) / da_c.TLARs.nengine;    % OEI
WoWTO(4)      = 1;                  % MTOM
CL = [CL,CLmaxes_CR(:)];           % Cruise config.
% 5 - approach
CGD_norm(5)   = 0.021;              % Req. Climb Gradient
VoVmin(5)     = 1.5;                % Req V in terms of V/V_stall_L
dCD0_flaps(5) = da_c.TLARs.dCD0_f_App;         % CD0 increment due to flaps in APPROACH config
dCD0_LG(5)    = da_c.TLARs.dCD0_lgs;           % Land. Gear Extended 
e_vet(5)      = da_c.TLARs.e + da_c.TLARs.de_LND;
eng_w(5)      = (da_c.TLARs.nengine-1) / da_c.TLARs.nengine;    % OEI
WoWTO(5)      = da_c.MLndoMTo;         % Max. Landing Weight
CL = [CL,( CLmaxes_LND(:)+CLmaxes_TO(:) )*0.5]; % Approach as an average between LND and T/O 
% 6 - Balked Lading AEI
CGD_norm(6)   = 0.032;              % Req. Climb Gradient
VoVmin(6)     = 1.3;                % Req V in terms of V/V_stall_L
dCD0_flaps(6) = da_c.TLARs.dCD0_f_LND;         % CD0 increment due to flaps in LND config
dCD0_LG(6)    = da_c.TLARs.dCD0_lgs;           % Land. Gear Extended
e_vet(6)      = da_c.TLARs.e+da_c.TLARs.de_LND;
eng_w(6)      = 1;                  % AEI
WoWTO(6)      = da_c.MLndoMTo;         % Max. Landing Weight
CL = [CL,CLmaxes_LND(:)];          % LND 

COND = {'First Segment CL$_{max,TO}$ = ','Transition CL$_{max,TO}$ = ','Second Segment CL$_{max,TO}$ = ','En-Route CL$_{max}$ = ',...
    'Approach CL$_{max,LND}$ = ','Balked CL$_{max,app}$ = '};
hold( ax1,'on'); hold( ax2,'on');
for i = 1:nCLs
    %CL = CLmaxes(i)/(1.2^2);
    
    for j = 1:6
        K = 1/( pi*da_c.ARw*e_vet(j));
        CD = CD0 + dCD0_flaps(j) + dCD0_LG(j) + K*( CL(i,j)/VoVmin(j)^2 )^2;
        E = CL(i,j)/VoVmin(j)^2/CD;         % Efficiency
        ToW = CGD_norm(j) + 1/E;
        ToW = WoWTO(j) * ToW*da_c.TLARs.T0oTmc * TisaoT50 / eng_w(j);
        % {Kg/m^2]
        lin(i,j) = plot( ax1,[0,1000],ToW*[1,1] );
        lin(i,j).LineStyle = '-.'; lin(i,j).LineWidth = 2;
        lin(i,j).Color = colors(i,:);
        lin(i,j).DisplayName = [ COND{j},num2str(CL(i,j)) ];

        % {lb/ft^2]
        lin(i,j) = plot( ax2,[0,1000*2.204623/(3.28084^2)],ToW*[1,1] );
        lin(i,j).LineStyle = '-.'; lin(i,j).LineWidth = 2;
        lin(i,j).Color = colors(i,:);
        lin(i,j).DisplayName = [ COND{j},num2str(CL(i,j)) ];

    end

end

%% Additional Conditions
% Ceiling Requirement
if nargin > 8 && ~isempty( Roc_req )
    Emax = sqrt( 0.25*pi*da_c.ARw*da_c.TLARs.e/CD0 );    % Max Efficiency
    CLe  = sqrt( pi*da_c.ARw*da_c.TLARs.e*CD0 );         % CL at max Efficiency
    nRoCs = length(Roc_req(:,1));
    % Cruise Ceiling check
    for k = 1:nRoCs
        % Max RoC at Cruise ceiling must be 300 ft/min
        RoC         = Roc_req(k,1)*0.3048/60;           % Required RoC in [m/s]
        WoS         = linspace(0,1000,100);             % WoS in [kg/m^2]
        [T,a,P,rho] = atmosisa( Roc_req(k,2) );
        ToT0        = 1/( 0.75*(rho/1.225)*1 );         % Altitude Effect on Thrust at max Admission
        WoWTO       = [da_c.MCroMTo(1),da_c.MCroMTo(2),0.5*( da_c.MCroMTo(1)+da_c.MCroMTo(2) )];
        WEIGHT      = {'Initial','Final','Avg'};

        i = nCLs+k;
        for j = 1:3
            ToW  = WoWTO(j)/ToT0.*( RoC./(sqrt( 2/(rho*CLe) ).*sqrt(WoS*9.81).*sqrt(WoWTO(j)) ) + 1/Emax );
            % {Kg/m^2]
            lin(i,j) = plot( ax1,WoS,ToW);
            lin(i,j).LineStyle = '-.'; lin(i,j).LineWidth = 2;
            lin(i,j).DisplayName = [ 'ROC of ',num2str(RoC),' m/s at ',num2str(da_c.TLARs.cruise.h),'m; CL$_E$: ',num2str(CLe),' and ',WEIGHT{j},' Cruise Weight' ];
            lin(i,j).Color = colors(i,:);
            % [lb/ft^2]
            lin(i,j) = plot( ax2,WoS*2.204623/(3.28084^2),ToW);
            lin(i,j).LineStyle = '-.'; lin(i,j).LineWidth = 2;
            lin(i,j).DisplayName = [ 'Cruise Ceiling at ',num2str(da_c.TLARs.cruise.h),'m and ',num2str(CLe),' and ',WEIGHT{j} ];
            lin(i,j).Color = colors(i,:);
        end
    end
end
end