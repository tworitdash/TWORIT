function [] = plot_FF(plot_info)
    
    if plot_info.same_diff == 0
        figure; plot(plot_info.x, plot_info.y, 'DisplayName', plot_info.lgd, 'LineWidth', 2); 
        grid on;
    else
        hold on; plot(plot_info.x, plot_info.y, 'DisplayName', plot_info.lgd, 'LineWidth', 2); 
        grid on;
    end
    
    xlabel(plot_info.xl, 'FontSize', 12, 'FontWeight', 'bold');
    ylabel(plot_info.yl, 'FontSize', 12, 'FontWeight', 'bold');
    legend;
end