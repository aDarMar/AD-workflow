function ax_pl = sizing_plot_cfr( AirData,nAir,ax_pl )
%sizing_plot_cfr Function that places the given airplanes on the sizing
%plot
%   Detailed explanation goes here

%cf_ax = axes( 'Parent',fig_ri ); hold( cf_ax,'on' );

j = 1;
for i=1:nAir
    try
        if  ~isnan( AirData(i).ToW )
            lines_plot(j)             = plot( ax_pl,AirData(i).WoS, AirData(i).ToW );
            lines_plot(j).LineStyle   = 'none'; 
            lines_plot(j).Marker      = AirData(i).Mark;
            lines_plot(j).MarkerSize  = 4; lines_plot(j).LineWidth = 2;
            lines_plot(j).DisplayName = [ AirData(i).Family,' ',AirData(i).Name,'-',AirData(i).Code];
            j = j +1;
        end
    catch
        
    end
end

end

