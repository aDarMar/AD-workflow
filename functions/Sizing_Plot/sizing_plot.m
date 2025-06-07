function [ siz_ax,rem_ch,Roc_s ] = sizing_plot( da_c,iS,CLmax_TO_vett,...
    CLmax_LND_vett,CLmax_CR_vett,sigma,TisaoT50,V_cr_vett,h_cr_vett,phi_v, fig_ri,fig_aux,Roc_s,rem_ch )
%sizing_plot Summary of this function goes here
%   Detailed explanation goes here
%INPUT
%   rem_ch: vector with indices of choices made for plotting
%   Roc_s: vector whose rows are RoCs desidered ( in [ft7min at given height [m])
%OUTPUT
%   s_idx: vector containing the indices of the choices made 
%% Input Data
CD0 = da_c.SizHis(iS).CD0; dCD0_wave = da_c.dCD0_wave; 
WLNDoWTO = da_c.MLndoMTo; WcroWTO = da_c.MCroMTo;
%% Graphics
% Initializzation
delete( fig_aux.Children )
aus_ax_1 = subplot(2,1,1,'Parent',fig_aux); aus_ax_2 = subplot(2,1,2,'Parent',fig_aux);
% Plotting
%[ aus_ax_1,aus_ax_2 ] = sizing_plot_TO( aus_ax_1,aus_ax_2,da_c.TLARs.TO.fieldmax,CLmax_TO_vett,sigma );
% [ aus_ax_1,aus_ax_2 ] = sizing_plot_LND( aus_ax_1,aus_ax_2,da_c.TLARs.LND.SGmax,CLmax_LND_vett,sigma,WLNDoWTO );
% [ aus_ax_1,aus_ax_2 ] = sizing_plot_Climb( aus_ax_1,aus_ax_2,CLmax_TO_vett,CLmax_LND_vett,CLmax_CR_vett,...
%                             da_c,CD0,TisaoT50 );
% [ aus_ax_1,aus_ax_2 ]   = sizing_plot_Cruise( aus_ax_1,aus_ax_2,CD0,dCD0_wave,V_cr_vett,h_cr_vett,WcroWTO,da_c.TLARs.e,da_c.ARw,phi_v );

%% Plot Choices
% Take-Off
nTO = length(CLmax_TO_vett);
n_inp = 14; % excludes the last input rem_ch
if nargin < n_inp
    rem_ch = nan(1,2); % Initializes rem_ch
    disp( '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    disp( '%%%%%%%%%%%%%%%%%%%%%% TAKE-OFF %%%%%%%%%%%%%%%%%%%%%%%%')
    disp(' Choose CL max at TO to display' );
    disp( CLmax_TO_vett );
    scelte = scelta_fun(nTO);
    rem_ch = idx_fun( rem_ch,scelte,'T/O' );
end
vet_idx = 1*ones(nTO,1);

% Landing
nTO = length(CLmax_LND_vett);
if nargin < n_inp
    disp( '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    disp( '%%%%%%%%%%%%%%%%%%%%%% LANDING %%%%%%%%%%%%%%%%%%%%%%%%')
    disp(' Choose CL max at LANDING to display' );
    disp( CLmax_LND_vett );
    scelte = scelta_fun(nTO);
    rem_ch = idx_fun( rem_ch,scelte,'LND' );
end
vet_idx = [vet_idx;2*ones(nTO,1)];

% CLIMB
nTO = length(CLmax_LND_vett)*length(CLmax_CR_vett)*length(CLmax_TO_vett); % max possible iterations
CLi = {'CL@T/O','CL@Cr','CL@LND'};
temp = [CLmax_TO_vett;CLmax_CR_vett;CLmax_LND_vett]; ct = 1;
if nargin < n_inp
    disp( '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    disp( '%%%%%%%%%%%%%%%%%%%%%%%% CLIMB %%%%%%%%%%%%%%%%%%%%%%%%%%')
    disp(' Choose the combination of CL maxes at each Climb phase: ' );
    while ct < nTO+1
        disp(['---- Condtion: ',num2str(ct),' ----'] )
        k = 1; tpm = nan(4,1); % Vector that stores all the input coices
        scelte = 1*nTO;
        while k < 4 && scelte > 0
            tpm(k) = scelte;
            disp( [CLi{k},'>> ',num2str(temp(k,:)) ] );
            scelte = scelta_fun(1);
            k = k +1;
        end
        tpm(end) = scelte;
        tpm = tpm(2:4); % excluding the first element that is nTO by definition
        if length( tpm( tpm>0 ) ) < 3
            % This means that an invalid number has been inserted to
            % terminate the sequence
            break
        else
            rem_ch = idx_fun( rem_ch,tpm,'Climb' );
            ct = ct + 1;
        end
    end
end
climb_idxs = rem_ch( rem_ch(:,1) == 3,2 ); % Indices associated to Climb
climb_idxs = reshape( climb_idxs,3, length( climb_idxs )/3 )'; % reshapes climb_idxs as a n_connd x 3 matrix
% In this way climb_idxs will be:
%   [idx_CL_TO,idx_CL_cruise,idx_CL_LAND]
vet_idx = [ vet_idx;3*ones( 6*length(climb_idxs(:,1)),1 ) ];
% CLIMB PT. 2
if nargin < n_inp-1
    disp( '%%%%%%%%%%%%%%%%%%%%%%%% CLIMB %%%%%%%%%%%%%%%%%%%%%%%%%%')
    disp(' Choose the ROC required at a fixed height: ' );
    j = 1; tmp = 1; Roc_s = [];
    while tmp > 0
        disp('ROC [ft/min]')
        tmp = input('>>');
        if tmp > -0.01
            Roc_s(j,1)  = tmp;
            disp('h [ft]')
            tmp = input('>>');
            if tmp > 0
                Roc_s(j,2) = tmp*0.3048;
                j = j + 1;
            else
                Roc_s = Roc_s(1:end-1,:);
                break
            end
        else
            break
        end
    end
end
if ~isempty( Roc_s )
    vet_idx = [ vet_idx;4*ones( 3*length(Roc_s(:,1)),1 ) ];
end

% CRUISE
nTO = length(h_cr_vett);
if nargin < n_inp
    disp( '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    disp( '%%%%%%%%%%%%%%%%%%%%%% CRUISE %%%%%%%%%%%%%%%%%%%%%%%%')
    disp(' Choose the h-V-phi combination in CRUISE to display' );
    disp( 1:nTO )
    disp( [' Altitude [m]   >> ',num2str(h_cr_vett)] );
    disp( [' Speed    [m/s] >> ',num2str(V_cr_vett)] );
    disp( [' Admission [%]  >> ',num2str(phi_v)] );
    scelte = scelta_fun(nTO);
    rem_ch = idx_fun( rem_ch,scelte,'Cruise' );
    rem_ch = rem_ch(2:end,:); % Removes the first row as it is a nan row
end
vet_idx = [vet_idx;rem_ch(end,1)*ones(nTO*3,1)];

disp( '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

%% Plot

[ aus_ax_1,aus_ax_2 ] = sizing_plot_TO( aus_ax_1,aus_ax_2,da_c.TLARs.TO.fieldmax,CLmax_TO_vett,sigma );
[ aus_ax_1,aus_ax_2 ] = sizing_plot_LND( aus_ax_1,aus_ax_2,da_c.TLARs.LND.SGmax,CLmax_LND_vett,sigma,WLNDoWTO );
[ aus_ax_1,aus_ax_2 ] = sizing_plot_Climb( aus_ax_1,aus_ax_2,...
    CLmax_TO_vett( climb_idxs(:,1) ),CLmax_LND_vett( climb_idxs(:,1) ),CLmax_CR_vett( climb_idxs(:,1) ),...
                            da_c,CD0,TisaoT50,Roc_s );
[ aus_ax_1,aus_ax_2 ] = sizing_plot_Cruise( aus_ax_1,aus_ax_2,CD0,dCD0_wave,V_cr_vett,h_cr_vett,WcroWTO,da_c.TLARs.e,da_c.ARw,phi_v );

delete( fig_ri.Children );
siz_ax = copyobj( aus_ax_1,fig_ri ); lin_pl = siz_ax.Children;



siz_ax.XLim = [0,1000]; siz_ax.YLim = [0,1]; legend( siz_ax,'Interpreter','Latex' );
siz_ax.Position = [0.1300 0.1100 0.7750 0.8150];
chs_UI( vet_idx ,lin_pl,rem_ch )

% nCH = length( lin_pl ); j = 1;
% for i = 1:nCH
%     if lin_pl(i).Visible == 1
%         rem_vs(j) = i;
%         j = j + 1;
%     end
% end

end

function chs = idx_fun( chs,scelte,cond )
%idx_fun: function that returns a struct with indices of chosen conditions
%to plot
%INPUT
%   chs: struct of chosen conditions to plot 
%       chs.idxs = index chosen to plot
%       chs.cond = flight condition associated
%       scelte: vector of idices given in input
switch cond
    case 'T/O'
        cond = 1;
    case 'LND'
        cond = 2;
    case 'Climb'
        cond = 3;
    case 'Cruise'
        cond = 4;
    otherwise
        error('Condition given not recognised')
end

len = length( scelte );
len2 = length( chs(:,1) );
for i = 1:len
    chs( len2+i,2 ) = scelte(i);
    chs( len2+i,1 ) = cond;
end

end

function scelta = scelta_fun(nTO)
    scelta = nan(nTO,1);
    disp('0 to end');
    i = 1; temp = 1;
    while i<nTO+1 && temp>0
        scelta(i) = temp;
        temp      = input('>>');
        i = i+1;
    end
    
    if nTO > 1
        scelta(i) = temp;
        scelta    = scelta(2:end-1);
    else
        scelta = temp;
    end
end

function chs_UI(vet_idx,lin,rem_chs)
nCH = length(vet_idx); %chs = 0.1; 
if nargin == 3
    n_cond = nan(5,1); sum = 0; n = 1; n_inp = length( rem_chs(:,1) )+1;
    k = 1;
    while k < 5
        n_cond(k)       = length ( vet_idx(vet_idx == k) );
        compl_idx_range = 1:n_cond(k);
        hide_idxs       = setdiff( compl_idx_range,rem_chs( rem_chs(:,1) == k,2) );
        for n = hide_idxs
            m = n;
            m = m + sum;
            m = nCH - m + 1; % lin is ordered from last plotted to first
            if  lin(m).Visible == 1
                lin(m).Visible = 0;
                lin(m).Annotation.LegendInformation.IconDisplayStyle = 'off';  % Rimuove dalla legenda
            % else
            %     lin(chs).Visible = 0;
            %     lin(chs).Annotation.LegendInformation.IconDisplayStyle = 'off';  % Rimuove dalla legenda
            end
        end
        sum = sum + n_cond(k);
        k   = k + 1;
    end

end
plot_UI(vet_idx,lin)
disp('Show what')
chs = input('>>');
chs = nCH - chs + 1; % lin is ordered from last plotted to first
while chs > 0 && chs < nCH+1
    if  lin(chs).Visible == 1
        lin(chs).Visible = 0;
        lin(chs).Annotation.LegendInformation.IconDisplayStyle = 'off';  % Rimuove dalla legenda
    else
        lin(chs).Visible = 1;
        lin(chs).Annotation.LegendInformation.IconDisplayStyle = 'on';  % Rimuove dalla legenda
    end
    plot_UI(vet_idx,lin)
    disp('Show what')
    chs = input('>>');
    chs = nCH - chs + 1; % lin is ordered from last plotted to first
end

end

function plot_UI(vet_idx,lin)
    k = 1; j = 1; %l = 1; 
    nCH = length(vet_idx);
    COND = {'Take-Off','Landing','Climb','Cruise'}; n_cond = length( COND );
    while k < n_cond+1
        disp( ['---------',COND{k},'---------'] )
        while vet_idx(j) == k && j < nCH
            if lin(nCH - j +1).Visible == 1
                disp([num2str(j),' : S |',lin(nCH - j +1).DisplayName])
            else
                disp([num2str(j),' : NS|',lin(nCH - j +1).DisplayName])
            end
            % if vet_chs(l,1) == vet_idx(j)
            %     
            %     l = l + 1;
            % else
            %     disp([num2str(j),' : NS|',lin(j).DisplayName])
            % end
            j = j + 1;
        end
    k = k + 1;
    end
    j = nCH;
    if lin(nCH - j +1).Visible == 1
        disp([num2str(j),' : S |',lin(nCH - j +1).DisplayName])
    else
        disp([num2str(j),' : NS|',lin(nCH - j +1).DisplayName])
    end
end