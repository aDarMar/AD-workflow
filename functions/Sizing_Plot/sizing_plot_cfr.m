function [f,LEG] = sizing_plot_cfr(AirData,nAir,fig,lines_plot,LEGin,siz_pt_line)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

figure(fig.Number)
if nargin == 6
    lines_plot = [lines_plot;siz_pt_line];
    LEGin = {LEGin{1:end},'Sizing Point'};
end

n_tmp = length(lines_plot);
j = 1;
for i=1:nAir
    try
        if  ~isnan( AirData(i).ToW )
            lines_plot(n_tmp+j) = plot( AirData(i).WoS, AirData(i).ToW );
            lines_plot(n_tmp+j).LineStyle = 'none'; lines_plot(n_tmp+j).Marker = AirData(i).Mark;
            lines_plot(n_tmp+j).MarkerSize = 4; lines_plot(n_tmp+j).LineWidth = 2;
            LEG{j} = [ AirData(i).Family,' ',AirData(i).Name,'-',AirData(i).Code] ;
            j = j +1;
        end
    catch
        
    end
end

n_tmp = length(lines_plot);
ax_fig = axes('Parent', fig); axis([0,1000,0,1]); hold on;
grid minor; xlabel('W$_{TO}$ / S [Kg/m$^2$]','Interpreter','latex','FontSize',16);
ylabel('T$_{TO}$ / W$_{TO}$ [-]','Interpreter','latex','FontSize',16);
for i = 1:n_tmp
    copyobj(lines_plot(i),ax_fig);
end
legend( ax_fig, { LEGin{1:end},LEG{1:end} },'Interpreter','latex','FontSize',16 )

end

