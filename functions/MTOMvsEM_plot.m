function  MTOMvsEM_plot(aero_obj,nAero,lin_reg,lin_reg_cfr)
%MTOMvsEM_plot crea il grafico dell' MTOM in funzione dell' Empty Weight
%   Detailed explanation goes here

%GRF_m = [223, 255, 0; 255, 191, 0; 100, 149, 237; ]

figure()



i = 1;
subplot(1,2,1)
fig(1,i) = plot( aero_obj(i).EM,aero_obj(i).MTOM ); hold on
fig(1,i).LineStyle = 'none'; fig(1,i).Marker = aero_obj(i).Mark; fig(1,i).MarkerSize = 4; fig(1,i).LineWidth = 2;
grid minor; xlabel( 'MTOM [Kg]' ); ylabel( 'EM [Kg]' ); axis equal
subplot(1,2,2)
fig(2,i) = loglog( aero_obj(i).MTOM,aero_obj(i).EM ); hold on
fig(2,i).LineStyle = 'none'; fig(2,i).Marker = aero_obj(i).Mark; fig(2,i).MarkerSize = 4; fig(2,i).LineWidth = 2;
leg{1,i} = [ aero_obj(i).Family,' ',aero_obj(i).Name,'-',aero_obj(i).Code] ;
grid minor; xlabel( 'MTOM [Kg]' ); ylabel( 'EM [Kg]' )

for i = 2:nAero
    subplot(1,2,1)
    fig(1,i) = plot( aero_obj(i).EM,aero_obj(i).MTOM );
    fig(1,i).LineStyle = 'none'; fig(1,i).Marker = aero_obj(i).Mark; fig(1,i).MarkerSize = 4; fig(1,i).LineWidth = 2;
    col = rand(1,3); fig(1,i).MarkerEdgeColor = col;
    subplot(1,2,2)
    fig(2,i) = loglog( aero_obj(i).MTOM,aero_obj(i).EM );
    fig(2,i).LineStyle = 'none'; fig(2,i).Marker = aero_obj(i).Mark; fig(2,i).MarkerSize = 4; fig(2,i).LineWidth = 2;
    fig(2,i).MarkerEdgeColor = col;
    leg{1,i} = [ aero_obj(i).Family,' ',aero_obj(i).Name,'-',aero_obj(i).Code] ;
end

if nargin > 2
    Wmtom_reg  = (2:20)*1e4; Wmtom_reg = Wmtom_reg*2.2046;              % da [Kg] -> [lb]
    Wempty_reg = 10.^( ( log10(Wmtom_reg) - lin_reg(1) )./lin_reg(2) ); % in [lb]
    Wempty_reg = Wempty_reg/2.2046; Wmtom_reg = Wmtom_reg/2.2046;       % da [lb] -> [Kg]

    i = nAero+1;
    subplot(1,2,1)
    fig(1,i) = plot( Wempty_reg,Wmtom_reg );
    fig(1,i).LineStyle = '-'; fig(1,i).LineWidth = 2; %fig(1,i).Marker = aero_obj(i).Mark; fig(1,i).MarkerSize = 4;
    col = rand(1,3); fig(1,i).MarkerEdgeColor = col;
    subplot(1,2,2)
    fig(2,i) = plot( Wmtom_reg,Wempty_reg );
    fig(2,i).LineStyle = '-'; fig(2,i).LineWidth = 2; %fig(1,i).Marker = aero_obj(i).Mark; fig(1,i).MarkerSize = 4;
    col = rand(1,3); fig(2,i).MarkerEdgeColor = col;
    leg{1,i} =  ['log_{10}(W_{MTOM}) = ',num2str( lin_reg(2) ),'log_{10}(W_{E})  ',num2str( lin_reg(1) )] ;
    
    if nargin == 4
        n_cases = length( lin_reg_cfr(:,1) );
        for j = 1:n_cases
            Wmtom_reg  = (2:20)*1e4; Wmtom_reg = Wmtom_reg*2.2046;              % da [Kg] -> [lb]
            Wempty_reg = 10.^( ( log10(Wmtom_reg) - lin_reg_cfr(j,1) )./lin_reg_cfr(j,2) ); % in [lb]
            Wempty_reg = Wempty_reg/2.2046; Wmtom_reg = Wmtom_reg/2.2046;       % da [lb] -> [Kg]
            i = nAero+1+j;
            subplot(1,2,1)
            fig(1,i) = plot( Wempty_reg,Wmtom_reg );
            fig(1,i).LineStyle = '-'; fig(1,i).LineWidth = 2; %fig(1,i).Marker = aero_obj(i).Mark; fig(1,i).MarkerSize = 4;
            col = rand(1,3); fig(1,i).MarkerEdgeColor = col;
            subplot(1,2,2)
            fig(2,i) = plot( Wmtom_reg,Wempty_reg );
            fig(2,i).LineStyle = '-'; fig(2,i).LineWidth = 2; %fig(1,i).Marker = aero_obj(i).Mark; fig(1,i).MarkerSize = 4;
            col = rand(1,3); fig(2,i).MarkerEdgeColor = col;
            leg{1,i} =  ['log_{10}(W_{MTOM}) = ',num2str( lin_reg_cfr(j,2) ),'log_{10}(W_{E})  ',num2str( lin_reg_cfr(j,1) )] ;
            
        end
    end
    
end

subplot(1,2,1)
legend( fig(1,:),leg{1,:} )
subplot(1,2,2)
legend( fig(2,:),leg{1,:} )

end

