function plotMaxMinRowNorm_degreeTotal(vRatio_MaxMin, myLimits, limits)
%
% % Inputs
%
% vRatio_MaxMin : 
% 
% myLimits : 
%
% limits :



myLowerLimit = myLimits(1);
myUpperLimit = myLimits(2);

lowerLimit = limits(1);
upperLimit = limits(2);

figure_name = sprintf([mfilename ' : ' 'Plot Max : Min Row Norm']);
figure('name',figure_name)
hold on
vec_x = myLowerLimit : 1 : myUpperLimit;
plot(vec_x, log10(vRatio_MaxMin),'-s')
vline(lowerLimit)
vline(upperLimit)
hold off

end
