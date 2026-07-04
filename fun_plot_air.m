function   fun_plot_air(Air,nAir)
%FUN_PLOT_AIR Summary of this function goes here
%   Detailed explanation goes here

fig = figure('Name','Geometrical Data from Reference Aircrafts');

ax1 = subplot(2,2,1,'Parent',fig); ax2 = subplot(2,2,2,'Parent',fig);
ax3 = subplot(2,2,3,'Parent',fig); 
hold( ax1,'on' ); hold( ax2,'on' ); hold( ax3,'on' ); 

for i = 1:nAir
    col = rand(1,3);
    j = 1; lin(i,j) = plot( ax1,Air(i).wing.bw, Air(i).wing.panels(1).TR );
    lin(i,j).LineStyle = 'none'; lin(i,j).Marker = Air(i).Mark; lin(i,j).MarkerSize = 4; lin(i,j).LineWidth = 2;
    lin(i,j).MarkerEdgeColor = col;
    j = 2; lin(i,j) = plot( ax1,Air(i).wing.bw, Air(i).wing.TR );
    lin(i,j).LineStyle = 'none'; lin(i,j).Marker = Air(i).Mark; lin(i,j).MarkerSize = 4; lin(i,j).LineWidth = 2;
    lin(i,j).MarkerEdgeColor = col;
    j = 3; lin(i,j) = plot( ax2,Air(i).wing.Sw, 2*Air(i).wing.panels(1).S/Air(i).wing.Sw );
    lin(i,j).LineStyle = 'none'; lin(i,j).Marker = Air(i).Mark; lin(i,j).MarkerSize = 4; lin(i,j).LineWidth = 2;
    lin(i,j).MarkerEdgeColor = col;
    j = 4; lin(i,j) = plot( ax3,0.5*Air(i).wing.bw, 2*Air(i).wing.panels(1).tip.yglob/Air(i).wing.bw );
    lin(i,j).LineStyle = 'none'; lin(i,j).Marker = Air(i).Mark; lin(i,j).MarkerSize = 4; lin(i,j).LineWidth = 2;
    lin(i,j).MarkerEdgeColor = col;
end

end

