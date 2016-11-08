
switch SETTINGS.PLOT_GRAPHS
    case 'y'
        
        title = sprintf('%s - Residuals',mfilename());
        xlabel('Iteration Number')
        ylabel('log_{10} Residual')
        figure('name',title)
        hold on
        plot(log10(condition),'-s')
        hold off
    case 'n'
end
