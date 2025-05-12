function  [fig,LEG] = sizing_plot_Climb( ...
    CLmaxes_TO,CLmaxes_LND,CLmaxes_CR,e,de_TO,de_LND,...
    CD0,dCD0_f_TO,dCD0_f_LND,dCD0_f_APP,dCD0_lgs,AR,...
    neng,WlandoMTOM,T0oTmc,TisaoT50 )
%sizing_plot_Climb Function for plotting climb requirements
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

nCLs = length(CLmaxes_TO); 
%CL = [CLmaxes_TO(:)';CLmaxes_LND(:)';CLmaxes_CR(:)'];
% Conditions:
% 1 - first segment of climb
CGD_norm(1)   = 0.012;              % Req. Climb Gradient
VoVmin(1)     = 1.2;                % Req V in terms of V/V_stall_TO
dCD0_flaps(1) = dCD0_f_TO;          % CD0 increment due to flaps in T/O config
dCD0_LG(1)    = 0;                  % Land. Gear Retracted
e_vet(1)          = e+de_TO;
eng_w(1)      = (neng-1) / neng;    % OEI
WoWTO(1)      = 1;                  % MTOM
CL = CLmaxes_TO(:);                 % T/O 
% 2 - transitio to climb
CGD_norm(2)   = 0;                  % Req. Climb Gradient
VoVmin(2)     = 1.1;                %Req V in terms of V/V_stall_TO
dCD0_flaps(2) = dCD0_f_TO;          % CD0 increment due to flaps in T/O config
dCD0_LG(2)    = dCD0_lgs;           % Land. Gear Extended
e_vet(2)      = e+de_TO;
eng_w(2)      = (neng-1) / neng;    % OEI
WoWTO(2)      = 1;                  % MTOM
CL = [CL,CLmaxes_TO(:)];            % T/O 
% 3 - second segment
CGD_norm(3)   = 0.024;              % Req. Climb Gradient
VoVmin(3)     = 1.2;                % Req V in terms of V/V_stall_L
dCD0_flaps(3) = dCD0_f_TO;          % CD0 increment due to flaps
dCD0_LG(3)    = 0;                  % Land. Gear Extended
e_vet(3)      = e+de_TO;
eng_w(3)      = (neng-1) / neng;    % OEI
WoWTO(3)      = 1;                  % MTOM
CL = [CL,CLmaxes_TO(:)]; % T/O 
% 4 - en-route Climb
CGD_norm(4)   = 0.012;              % Req. Climb Gradient
VoVmin(4)     = 1.25;               % Req V in terms of V/V_stall_L
dCD0_flaps(4) = 0;                  % CD0 increment due to flaps
dCD0_LG(4)    = 0;                  % Land. Gear Retracted
e_vet(4)      = e;
eng_w(4)      = (neng-1) / neng;    % OEI
WoWTO(4)      = 1;                  % MTOM
CL = [CL,CLmaxes_CR(:)];           % Cruise config.
% 5 - approach
CGD_norm(5)   = 0.021;              % Req. Climb Gradient
VoVmin(5)     = 1.5;                % Req V in terms of V/V_stall_L
dCD0_flaps(5) = dCD0_f_APP;         % CD0 increment due to flaps in APPROACH config
dCD0_LG(5)    = dCD0_lgs;           % Land. Gear Extended 
e_vet(5)      = e + de_LND;
eng_w(5)      = (neng-1) / neng;    % OEI
WoWTO(5)      = WlandoMTOM;         % Max. Landing Weight
CL = [CL,( CLmaxes_LND(:)+CLmaxes_TO(:) )*0.5]; % Approach as an average between LND and T/O 
% 6 - Balked Lading AEI
CGD_norm(6)   = 0.032;              % Req. Climb Gradient
VoVmin(6)     = 1.3;                % Req V in terms of V/V_stall_L
dCD0_flaps(6) = dCD0_f_LND;         % CD0 increment due to flaps in LND config
dCD0_LG(6)    = dCD0_lgs;           % Land. Gear Extended
e_vet(6)      = e+de_LND;
eng_w(6)      = 1;                  % AEI
WoWTO(6)      = WlandoMTOM;         % Max. Landing Weight
CL = [CL,CLmaxes_LND(:)];          % LND 

COND = {'First Segment CL$_{max,TO}$ = ','Transition CL$_{max,TO}$ = ','Second Segment CL$_{max,TO}$ = ','En-Route CL$_{max}$ = ',...
    'Approach CL$_{max,LND}$ = ','Balked CL$_{max,app}$ = '};
for i = 1:nCLs
    %CL = CLmaxes(i)/(1.2^2);
    
    for j = 1:6
        K = 1/( pi*AR*e_vet(j));
        CD = CD0 + dCD0_flaps(j) + dCD0_LG(j) + K*( CL(i,j)/VoVmin(j)^2 )^2;
        E = CL(i,j)/VoVmin(j)^2/CD;         % Efficiency
        ToW = CGD_norm(j) + 1/E;
        ToW = WoWTO(j) * ToW*T0oTmc * TisaoT50 / eng_w(j);
        
        subplot 212; hold on
        fig(i,j) = plot( [0,200],ToW*[1,1] );
        fig(i,j).LineStyle = '-.'; fig(i,j).LineWidth = 2; %fig(1,i).Marker = aero_obj(i).Mark; fig(1,i).MarkerSize = 4;
        col = rand(1,3); fig(i,j).MarkerEdgeColor = col; 
        %set(fig(i,j), 'Visible','off');
        subplot 211; hold on
        fig(i,j) = plot( [0,1000],ToW*[1,1] );
        fig(i,j).LineStyle = '-.'; fig(i,j).LineWidth = 2; %fig(1,i).Marker = aero_obj(i).Mark; fig(1,i).MarkerSize = 4;
        col = rand(1,3); fig(i,j).MarkerEdgeColor = col; 
        %set(fig(i,j), 'Visible','off');
        LEG{i,j} = [ COND{j},num2str(CL(i,j)) ];
        
    end

end
end