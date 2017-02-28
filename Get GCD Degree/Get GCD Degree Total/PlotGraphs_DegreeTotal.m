global SETTINGS
if(SETTINGS.PLOT_GRAPHS)
    
    
    plot_R_diagonals = false;
    plot_MinSingularValues = true;
    plot_QRMaxMinDiags = false;
    plot_QRMaxMinRowSum = false;
    
    x = lower_lim_comp:1:upper_lim_comp;
    
    if plot_R_diagonals == true
        figure_name = sprintf('%s - R Diagonals',mfilename);
        figure('name',figure_name)
        hold on
        title('Diagonal entries of R matrix where S_{k,k} = QR')
        xlabel('k')
        ylabel('log_{10} diagonal r_{i,i}')
        xlim([1,+inf])
        scatter(data(:,1),log10(data(:,2)))
        vline(lower_lim,'b','');
        vline(upper_lim,'b','');
        hold off
        MySave('TotalDegree_RDiagonals');
    end
    
    if plot_MinSingularValues == true
        figure_name = sprintf('%s - Minimum Singular Values',mfilename);
        figure('name',figure_name)
        titleString = sprintf(['Minimal singular values of each subresultant S_{k} \n']);
        title(titleString)
        xlabel('k: total degree')
        ylabel('log_{10} \sigma_{i}')
        hold on
        xlim([1,+inf])
        plot(x,log10(v_MinimumSingularValue),'-s');
        vline(lower_lim,'b','');
        vline(upper_lim,'b','');
        hold off
        MySave('TotalDegree_MinimumSingularValues');
    end
    
    if plot_QRMaxMinDiags == true
        % plot all the largest ratios for k = 1,...,min(m,n)
        figure_name = sprintf('%s - QR max:min diagonals',mfilename);
        figure('name',figure_name)
        hold on
        title('Plotting max:min diagonal entries of QR decomposition of S_{k}')
        vRatio_MaxMin_Diags_R1 = v_maxDiagR1 ./ v_minDiagR1;
        plot(x,log10(vRatio_MaxMin_Diags_R1),'-s');
        xlabel('k')
        ylabel('log_{10} max r_{k} / min r_{k}')
        xlim([1,+inf])
        vline(lower_lim,'b','');
        vline(upper_lim,'b','');
        hold off
        MySave('TotalDegree_RatioMaxMinDiagonals');
    end
    
    if plot_QRMaxMinRowSum == true
        figure_name = sprintf('%s - QR max:min Row Sum',mfilename);
        figure('name',figure_name)
        hold on
        title('Plotting max:min rowsums of QR decomposition of S_{k}')
        xlim([1,+inf])
        vRatio_MaxMin_RowNorm_R1 = v_maxRowNormR1 ./ v_minRowNormR1;
        plot(x,log10(vRatio_MaxMin_RowNorm_R1),'-s');
        xlabel('k')
        ylabel('log_{10} max r_{k} / min r_{k}')
        vline(lower_lim,'b','');
        vline(upper_lim,'b','');
        hold off
        MySave('TotalDegree_RatioMaxMinRowSum');
    end
    
    
end