global SETTINGS

% Get the name of the function which called this function.
[Stack,~] = dbstack();
calling_function = Stack(2).name;


if(SETTINGS.PLOT_GRAPHS)
    % Plot the minimum singular values of S_{k_{1},k_{2}} for all valid k1,k2
    % pairs.
    
    figure_name = sprintf('%s : Minimum Singular Values',calling_function);
    figure('name',figure_name)
    hold on
    xlim([0,+inf])
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
    
    figure()
    y = 0:1:min(m1,n1);
    x = 0:1:min(m2,n2);
    [X,Y] = meshgrid(x,y);
    F = log10(mat_MinimumSingularValues);
    surf(X,Y,F);
    hold off
    

    % Print the singular values
    try
        figure_name = sprintf('%s : Singular Values',calling_function);
        figure('name',figure_name)
        hold on
        xlim([1,+inf])
        plot(log10(vMinimumSingularValues),'-s')
        hold off
        MySave('DegreeRelative_SingularValues');
    catch
    end
    
end