function [] = plotMinimumSingularValues_degreeRelative(mat_MinimumSingularValues, limits_t1, limits_t2)
%
% % Inputs
%
% mat_MinimumSingularValues :
%
% limits_t1 :
% 
% limits_t2 :

lowerLimit_k1 = limits_t1(1);
upperLimit_k1 = limits_t1(2);
lowerLimit_k2 = limits_t2(1);
upperLimit_k2 = limits_t2(2);

figure_name = 'Minimum Singular Values';
figure('name', figure_name)

hold on

v_k1 = lowerLimit_k1 : 1 : upperLimit_k1;
v_k2 = lowerLimit_k2 : 1 : upperLimit_k2;

[X,Y] = meshgrid(v_k2,v_k1);

mesh(X,Y,log10(mat_MinimumSingularValues));

xlabel('k_{2}')
ylabel('k_{1}')
zlabel('log_{10} \sigma_{k1,k2}')

hold off

end