function [] = plotRowDiag_degreeTotal(arr_R1_diag, limits)

lowerLimit = limits(1);
upperLimit = limits(2);

nSubresultants = upperLimit - lowerLimit + 1;


figure_name = sprintf('Plot Row Diagonals of S');
figure('name',figure_name)

hold on
for i = 1:1:nSubresultants
   
    vR1_diags = arr_R1_diag{i};
    
    v_is = i.*ones(length(vR1_diags),1);
    
    plot(v_is,log10(vR1_diags),'*');
    
    
    
end
hold off

end