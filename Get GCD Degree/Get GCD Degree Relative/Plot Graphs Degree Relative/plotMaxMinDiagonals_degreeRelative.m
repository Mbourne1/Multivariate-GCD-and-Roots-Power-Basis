function [] = plotMaxMinDiagonals_degreeRelative(max_diags_R1, min_diags_R1, limits_t1, limits_t2)

lowerLimit_t1 = limits_t1(1);
upperLimit_t1 = limits_t1(2);
lowerLimit_t2 = limits_t2(1);
upperLimit_t2 = limits_t2(2);


ratio = (log10(min_diags_R1) ./ log10(max_diags_R1));
log10(ratio);

figure_name = '';
figure('name',figure_name)
hold on
x = lowerLimit_t2:1:upperLimit_t2;
y = lowerLimit_t1:1:upperLimit_t1;
[X, Y] = meshgrid(x,y);
surf(X,Y,(ratio));

hold off

end

