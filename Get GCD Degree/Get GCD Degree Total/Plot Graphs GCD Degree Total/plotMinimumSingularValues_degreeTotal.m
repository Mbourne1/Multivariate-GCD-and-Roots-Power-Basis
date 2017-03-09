function plotMinimumSingularValues_degreeTotal(vMinimumSingularValues, myLimits_t, limits_t)
%
% % Inputs
%
% vMinimumSingularValues :
%
% myLimits_t :
%
% limits_t :

myLowerLimit = myLimits_t(1);
myUpperLimit = myLimits_t(2);

lowerLimit = limits_t(1);
upperLimit = limits_t(2);

figure_name = sprintf([mfilename ' : ' 'Minimum Singular Values']);
figure('name',figure_name);
hold on

xlabel('k');
ylabel('log_{10} Singular Values');
x_vec = myLowerLimit : 1 : myUpperLimit;
plot(x_vec, log10(vMinimumSingularValues), '-s');
vline(lowerLimit)
vline(upperLimit)
hold off

end