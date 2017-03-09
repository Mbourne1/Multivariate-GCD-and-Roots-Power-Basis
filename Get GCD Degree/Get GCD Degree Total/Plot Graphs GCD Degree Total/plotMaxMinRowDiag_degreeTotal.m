function [] = plotMaxMinRowDiag_degreeTotal(vRatio_MaxMin_Diags_R1, myLimits, limits)
%
% % Inputs
%
% vRatio_MaxMin_Diags_R1
%
% myLimits : My limits on degree of GCD
%
% limits : Limits on degree of GCD


myLowerLimit = myLimits(1);
myUpperLimit = myLimits(2);

lowerLimit = limits(1);
upperLimit = limits(2);



figure_name = sprintf([mfilename ' : ' 'Plot Ratio Max over Min Diagonals of R1']);
figure('name',figure_name)
hold on
vec_x = myLowerLimit : 1 : myUpperLimit;
plot(vec_x, log10(vRatio_MaxMin_Diags_R1),'-s');
vline(lowerLimit);
vline(upperLimit);
hold off

end