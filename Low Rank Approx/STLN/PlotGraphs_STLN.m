
switch SETTINGS.PLOT_GRAPHS
    case 'y'
        
        title = sprintf('%s - Residuals',mfilename());
        figure('name',title)
        hold on
        plot(log10(condition),'-s')
        hold off
    case 'n'
end
