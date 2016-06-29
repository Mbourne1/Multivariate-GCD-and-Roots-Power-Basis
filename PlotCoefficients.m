function [] = PlotCoefficients(fxy,fww,name)
% PlotCoefficients(fxy,fww,name)
%
% Plot the coefficients of the polynomial f(x,y) and f(w,w)
%
% fxy : Matrix of coefficients of polynomial f(x,y)
%
% fww : Matrix of coefficients of polynomial f(w,w)
%
% name : Name of the function 'f'

global SETTINGS

switch SETTINGS.PLOT_GRAPHS
    case 'y'
        
        % Get Degree of polynomial f(x,y)
        [m1,m2] = GetDegree(fxy);
        
        % Get total number of coefficients in f(x,y)
        num_coeff_f = (m1 + 1) * (m2 + 1);
        
        % Get vector of coefficients of f(x,y)
        v_fxy = GetAsVector(fxy);
        
        % Get vector of coefficients of f(w,w)
        v_fww = GetAsVector(fww);
        
        label1 = sprintf('%s(x,y)',name);
        label2 = sprintf('%s(w,w)',name);
        
        title_label = sprintf('- Coefficients %s(x,y)',name);
        
        x_axis_f = (1:1:num_coeff_f);
        figure('name',strcat(mfilename(),title_label))
        hold on
        plot(x_axis_f, real(log10(v_fxy)), '-s','DisplayName',label1)
        plot(x_axis_f, real(log10(v_fww)),'-s','DisplayName',label2)
        
        ylabel('log_{10} |a_{i}|')
        xlabel('i')
        legend(gca,'show')
        hold off
        
    case 'n'
end
end