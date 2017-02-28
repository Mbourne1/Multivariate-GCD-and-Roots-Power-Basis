function [] = plotDiagonalsR1_degreeRelative(arr_DiagonalsR1, limits_t1, limits_t2)

lowerLimit_t1 = limits_t1(1);
upperLimit_t1 = limits_t1(2);
lowerLimit_t2 = limits_t2(1);
upperLimit_t2 = limits_t2(2);

figure_name = sprintf('');
figure('name',figure_name)
hold on

for i1 = lowerLimit_t1 : 1 : upperLimit_t1
    for i2 = lowerLimit_t2 : 1 : upperLimit_t2
    
        vec = abs(arr_DiagonalsR1{i1+1,i2+1});
        v_i1 = i1.* ones(length(vec));
        v_i2 = i2.* ones(length(vec));
        
        plot3(v_i1, v_i2, log10(vec), '*');
    end
end
grid on
hold off

end