function plotMaxMinRowNorm_degreeTotal(vRatio_MaxMin, limits)
%
% % Inputs
%
% vRatio_MaxMin
% 
% limits



figure_name = sprintf('Plot Max : Min Row Norm');
figure('name',figure_name)
hold on

plot(log10(vRatio_MaxMin),'-s')

hold off
