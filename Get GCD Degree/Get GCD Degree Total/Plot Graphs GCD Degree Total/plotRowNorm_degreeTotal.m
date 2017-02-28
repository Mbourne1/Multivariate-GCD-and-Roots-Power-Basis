function plotRowNorm_degreeTotal(arr_R1_RowNorm, limits)


lowerLimit = limits(1);
upperLimit = limits(2);
nSubresultants = upperLimit - lowerLimit + 1;


figure_name = sprintf('Row Norms');
figure('name',figure_name);
hold on

for i = 1:1:nSubresultants
    
    vR1RowNorm = arr_R1_RowNorm{i};
    
    % Produce a vector of i
    v_is = i.* ones(length(vR1RowNorm),1);
    
    plot(v_is,log10(vR1RowNorm),'*');
    
end
    
hold off


end