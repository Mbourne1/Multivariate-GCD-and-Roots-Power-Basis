function [] = plotMaxMinRowDiag_degreeTotal(vRatio_MaxMin_Diags_R1, limits)
%
%



figure_name = sprintf('Plot Ratio Max over Min Diagonals of R1');
figure('name',figure_name)
hold on
plot(log10(vRatio_MaxMin_Diags_R1),'-s');
hold off

end