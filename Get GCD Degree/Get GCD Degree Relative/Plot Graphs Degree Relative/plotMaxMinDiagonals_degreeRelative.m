function [] = plotMaxMinDiagonals_degreeRelative(max_diags_R1, min_diags_R1, myLimits_t1, myLimits_t2, limits_t1, limits_t2)
%
% % Inputs
%
% max_diags_R1 :
%
% min_diags_R1 :
%
% myLimits_t1 :
%
% myLimits_t2 :
%
% limits_t1 : 
%
% limits_t2 :

% Get upper and lower limits for graphing
myLowerLimit_t1 = myLimits_t1(1);
myUpperLimit_t1 = myLimits_t1(2);
myLowerLimit_t2 = myLimits_t2(1);
myUpperLimit_t2 = myLimits_t2(2);


ratio = (log10(min_diags_R1) ./ log10(max_diags_R1));



figure_name = sprintf([mfilename ' : ' 'MaxMin Row Diagonals' ]);
figure('name',figure_name)
hold on
x_vec = myLowerLimit_t2 : 1 : myUpperLimit_t2;
y_vec = myLowerLimit_t1 : 1 : myUpperLimit_t1;
[X, Y] = meshgrid(x_vec,y_vec);
surf(X,Y,(ratio));

hold off

end

