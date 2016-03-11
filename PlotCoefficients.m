switch PLOT_GRAPHS
    case 'y'
        % Get limit of the x axis for plotting coefficients of f and g
        x_axis_f = (1:1:num_coeff_f);
        x_axis_g = (1:1:num_coeff_g);
        
        % Plot the unprocessed and preprocessed coefficients of f(x,y) and
        % f(w,w)
        figure('name','Coefficients f(x,y)')
        hold on
        plot(x_axis_f, v_fxy_n, '-s','DisplayName','f(x,y)')
        plot(x_axis_f, v_fww_n,'-s','DisplayName','f(w,w)')
        legend(gca,'show')
        hold off
        
        % Plot the unprocessed and preprocessed coefficients of g(x,y) and
        % g(w,w)
        figure('name','Coefficients g(x,y)')
        hold on
        plot(x_axis_g, coeff_gxy_n,'-s','DisplayName','g(x,y) without processing')
        plot(x_axis_g, alpha.*coeff_gww_n,'-s','DisplayName','g(\omega_{1},\omega_{2})')
        legend(gca,'show')
        hold off
        
    case 'n'
    otherwise
        error('error - bool_plotgraphs is either y or n')
end