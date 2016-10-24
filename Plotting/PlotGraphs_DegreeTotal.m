global SETTINGS
switch SETTINGS.PLOT_GRAPHS
    case 'y'
        x = lower_lim_comp:1:upper_lim_comp;
        
        figure_name = sprintf('%s - R Diagonals',mfilename);
        figure('name',figure_name)
        hold on
        title('Diagonal entries of R matrix where S_{k,k} = QR')
        xlabel('k')
        ylabel('log_{10} diagonal r_{i,i}')
        scatter(data(:,1),log10(data(:,2)))
        vline(lower_lim,'b','');
        vline(upper_lim,'b','');
        hold off
        MySave('TotalDegree_RDiagonals');
        
        
        figure_name = sprintf('%s - Minimum Singular Values',mfilename);
        figure('name',figure_name)
        titleString = sprintf(['Minimal singular values of each subresultant S_{k} \n']);
        title(titleString)
        xlabel('k: total degree')
        ylabel('log_{10} \sigma_{i}')
        hold on
        plot(x,log10(v_MinimumSingularValue),'-s');
        vline(lower_lim,'b','');
        vline(upper_lim,'b','');
        hold off
        MySave('TotalDegree_MinimumSingularValues');
        
        % plot all the largest ratios for k = 1,...,min(m,n)
        figure_name = sprintf('%s - QR max:min diagonals',mfilename);
        figure('name',figure_name)
        hold on
        title('Plotting max:min diagonal entries of QR decomposition of S_{k}')
        vRatio_MaxMin_Diags_R1 = v_maxDiagR1 ./ v_minDiagR1;
        plot(x,log10(vRatio_MaxMin_Diags_R1),'-s');
        xlabel('k')
        ylabel('log_{10} max r_{k} / min r_{k}')
        vline(lower_lim,'b','');
        vline(upper_lim,'b','');
        hold off
        MySave('TotalDegree_RatioMaxMinDiagonals');
        
        figure_name = sprintf('%s - QR max:min Row Sum',mfilename);
        figure('name',figure_name)
        hold on
        title('Plotting max:min rowsums of QR decomposition of S_{k}')
        vRatio_MaxMin_RowNorm_R1 = v_maxRowNormR1 ./ v_minRowNormR1;
        plot(x,log10(vRatio_MaxMin_RowNorm_R1),'-s');
        xlabel('k')
        ylabel('log_{10} max r_{k} / min r_{k}')
        vline(lower_lim,'b','');
        vline(upper_lim,'b','');
        hold off
        MySave('TotalDegree_RatioMaxMinRowSum');
    case 'n'
    otherwise
        error('plot_graphs is either y or n')
end