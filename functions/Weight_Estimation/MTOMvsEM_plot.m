function  fig = MTOMvsEM_plot(aero_obj,nAero,lin_reg,lin_reg_cfr)
%MTOMvsEM_plot crea il grafico dell' MTOM in funzione dell' Empty Weight
%   Detailed explanation goes here

%GRF_m = [223, 255, 0; 255, 191, 0; 100, 149, 237; ]

fig = figure('Name','Weight Estimation');
ax_1 = subplot( 1,3,1,'Parent',fig ); ax_2 = subplot( 1,3,2,'Parent',fig );
i = 1;
lin(1,i) = plot( ax_1,aero_obj(i).EM,aero_obj(i).MTOM ); hold( ax_1, 'on' );
lin(1,i).LineStyle = 'none'; lin(1,i).Marker = aero_obj(i).Mark; lin(1,i).MarkerSize = 4; lin(1,i).LineWidth = 2;
ax_1.XGrid = 'on'; ax_1.YGrid = 'on'; axis(ax_1,'equal');
ax_1.XMinorGrid = 'on'; ax_1.YMinorGrid = 'on'; 
xlabel( 'MTOM [Kg]' ); ylabel( 'EM [Kg]' );

lin(2,i) = loglog( ax_2,aero_obj(i).MTOM,aero_obj(i).EM ); hold( ax_2, 'on' );
lin(2,i).LineStyle = 'none'; lin(2,i).Marker = aero_obj(i).Mark; lin(2,i).MarkerSize = 4; lin(2,i).LineWidth = 2;
leg{1,i} = [ aero_obj(i).Family,' ',aero_obj(i).Name,'-',aero_obj(i).Code] ;
ax_2.XGrid = 'on'; ax_2.YGrid = 'on'; %axis(ax_2,'equal');
ax_2.XMinorGrid = 'on'; ax_2.YMinorGrid = 'on'; 
xlabel( 'MTOM [Kg]' ); ylabel( 'EM [Kg]' )

max_MTOM = aero_obj(i).MTOM; min_MTOM = aero_obj(i).MTOM;
for i = 2:nAero
    % Find Maximum and Minimum MTOM for sizing the axes
    if max_MTOM < aero_obj(i).MTOM
       max_MTOM = aero_obj(i).MTOM; 
    end
    if min_MTOM > aero_obj(i).MTOM
        min_MTOM = aero_obj(i).MTOM;
    end
    % Plotting the Aircrafts
    % Linear Scale
    lin(1,i) = plot( ax_1,aero_obj(i).EM,aero_obj(i).MTOM );
    lin(1,i).LineStyle = 'none'; lin(1,i).Marker = aero_obj(i).Mark; lin(1,i).MarkerSize = 4; lin(1,i).LineWidth = 2;
    col = rand(1,3); lin(1,i).MarkerEdgeColor = col;
    % Logarithmic Scale
    lin(2,i) = loglog( ax_2,aero_obj(i).MTOM,aero_obj(i).EM );
    lin(2,i).LineStyle = 'none'; lin(2,i).Marker = aero_obj(i).Mark; lin(2,i).MarkerSize = 4; lin(2,i).LineWidth = 2;
    lin(2,i).MarkerEdgeColor = col;
    leg{1,i} = [ aero_obj(i).Family,' ',aero_obj(i).Name,'-',aero_obj(i).Code] ;
end

if nargin > 2
    % Plots the linear regression made with the input aircrafts
    Wmtom_reg  = linspace(min_MTOM,max_MTOM,10); Wmtom_reg = Wmtom_reg*2.2046;  % from [Kg] -> [lb]
    Wempty_reg = 10.^( ( log10(Wmtom_reg) - lin_reg(1) )./lin_reg(2) );         % in [lb]
    Wempty_reg = Wempty_reg/2.2046; Wmtom_reg = Wmtom_reg/2.2046;               % from [lb] -> [Kg]

    i = nAero+1;
    lin(1,i) = plot( ax_1,Wempty_reg,Wmtom_reg );
    lin(1,i).LineStyle = '-'; lin(1,i).LineWidth = 2; %lin(1,i).Marker = aero_obj(i).Mark; lin(1,i).MarkerSize = 4;
    col = rand(1,3); lin(1,i).MarkerEdgeColor = col;

    lin(2,i) = plot( ax_2,Wmtom_reg,Wempty_reg );
    lin(2,i).LineStyle = '-'; lin(2,i).LineWidth = 2; %lin(1,i).Marker = aero_obj(i).Mark; lin(1,i).MarkerSize = 4;
    col = rand(1,3); lin(2,i).MarkerEdgeColor = col;
    % Legend
    if lin_reg(1)>0
        leg{1,i} =  ['log_{10}(W_{MTOM}) = ',num2str( lin_reg(2) ),'log_{10}(W_{E})  + ',num2str( lin_reg(1) )] ;
    else
        leg{1,i} =  ['log_{10}(W_{MTOM}) = ',num2str( lin_reg(2) ),'log_{10}(W_{E})  - ',num2str( -lin_reg(1) )] ;
    end
    if nargin == 4
        n_cases = length( lin_reg_cfr(:,1) );
        for j = 1:n_cases
            Wmtom_reg  = linspace(min_MTOM,max_MTOM,10); Wmtom_reg = Wmtom_reg*2.2046;              % da [Kg] -> [lb]
            Wempty_reg = 10.^( ( log10(Wmtom_reg) - lin_reg_cfr(j,1) )./lin_reg_cfr(j,2) ); % in [lb]
            Wempty_reg = Wempty_reg/2.2046; Wmtom_reg = Wmtom_reg/2.2046;       % da [lb] -> [Kg]
            i = nAero+1+j;
            lin(1,i) = plot( ax_1,Wempty_reg,Wmtom_reg );
            lin(1,i).LineStyle = '-'; lin(1,i).LineWidth = 2; %lin(1,i).Marker = aero_obj(i).Mark; lin(1,i).MarkerSize = 4;
            col = rand(1,3); lin(1,i).MarkerEdgeColor = col;
            lin(2,i) = plot( ax_2,Wmtom_reg,Wempty_reg );
            lin(2,i).LineStyle = '-'; lin(2,i).LineWidth = 2; %lin(1,i).Marker = aero_obj(i).Mark; lin(1,i).MarkerSize = 4;
            col = rand(1,3); lin(2,i).MarkerEdgeColor = col;
            % Legend
            if lin_reg_cfr(j,1) < 0
                leg{1,i} =  ['log_{10}(W_{MTOM}) = ',num2str( lin_reg_cfr(j,2) ),'log_{10}(W_{E}) - ',num2str( -lin_reg_cfr(j,1) )] ;
            else
                leg{1,i} =  ['log_{10}(W_{MTOM}) = ',num2str( lin_reg_cfr(j,2) ),'log_{10}(W_{E}) + ',num2str( lin_reg_cfr(j,1) )] ;
            end
        end
    end
    
end

legend( ax_1,lin(1,:),leg{1,:} )
legend( ax_2,lin(2,:),leg{1,:} )

end

