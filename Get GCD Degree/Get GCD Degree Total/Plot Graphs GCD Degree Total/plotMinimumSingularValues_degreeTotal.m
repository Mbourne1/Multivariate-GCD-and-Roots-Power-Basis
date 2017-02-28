function plotMinimumSingularValues_degreeTotal(vMinimumSingularValues, limits)


figure_name = sprintf('Minimum Singular Values');
figure('name',figure_name);
hold on
plot(log10(vMinimumSingularValues), '-s');
hold off

end