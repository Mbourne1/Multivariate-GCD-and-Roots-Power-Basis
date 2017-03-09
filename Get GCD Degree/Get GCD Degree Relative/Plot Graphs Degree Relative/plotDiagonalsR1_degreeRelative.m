function [] = plotDiagonalsR1_degreeRelative(arr_DiagonalsR1, myLimits_t1, myLimits_t2, limits_t1, limits_t2)
%
% % Inputs
%
% arr_Diagonals_R1 :
%
% myLimits_t1 :
%
% myLimits_t2 :
%
% limits_t1 :
%
% limits_t2 :


lowerLimit_t1 = myLimits_t1(1);
upperLimit_t1 = myLimits_t1(2);
lowerLimit_t2 = myLimits_t2(1);
upperLimit_t2 = myLimits_t2(2);


nSubresultants_t1 = upperLimit_t1 - lowerLimit_t1 + 1;
nSubresultants_t2 = upperLimit_t2 - lowerLimit_t2 + 1;

figure_name = sprintf([mfilename ' : ' 'Diagonals of R1']);
figure('name',figure_name)
hold on

for i1 = 1 : 1 : nSubresultants_t1
    for i2 = lowerLimit_t2 : 1 : nSubresultants_t2
    
        k1 = lowerLimit_t1 + (i1-1);
        k2 = lowerLimit_t2 + (i2-1);
        
        vec = abs(arr_DiagonalsR1{i1+1, i2+1});
        
        v_k1 = k1.* ones(length(vec));
        v_k2 = k2.* ones(length(vec));
        
        plot3(v_k1, v_k2, log10(vec), '*');
        
    end
end
grid on
hold off

end