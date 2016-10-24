global SETTINGS

% Get the name of the function which called this function.
[Stack,~] = dbstack();
calling_function = Stack(2).name;


switch SETTINGS.PLOT_GRAPHS
    case 'y'
        % Plot the minimum singular values of S_{k_{1},k_{2}} for all valid k1,k2
        % pairs.
        figure_name = sprintf('%s : Minimum Singular Values',calling_function);
        figure('name',figure_name)
        hold on
        vMinimumSingularValues_all_log = log10(vMinimumSingularValues_all);
        scatter3(k1k2Pairs(:,1),k1k2Pairs(:,2),vMinimumSingularValues_all_log)
        xlabel('k1')
        ylabel('k2')
        zlabel('log_{10} \sigma')
        az = -45;
        el = 45;
        view(az, el);
        hold off
        MySave('DegreeRelative_MinimumSingularValues');
        
        
        
        %%
        % Print the singular values
        figure_name = sprintf('%s : Singular Values',calling_function);
        figure('name',figure_name)
        hold on
        plot(log10(vMinimumSingularValues),'-s')
        hold off
        MySave('DegreeRelative_SingularValues');
        
    case 'n'
    otherwise
        error('plot_graphs is either y or n')
end