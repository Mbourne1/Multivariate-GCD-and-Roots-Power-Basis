function [] = plotRowDiag_degreeTotal(arr_R1_diag, myLimits, limits)
%
% % Inputs
%
% arr_R1_diag :
%
% myLimits :
%
% limits

myLowerLimit = myLimits(1);
myUpperLimit = myLimits(2);

lowerLimit = limits(1);
upperLimit = limits(2);

nSubresultants = myUpperLimit - myLowerLimit + 1;

figure_name = sprintf([mfilename ' : ' 'Plot Row Diagonals of S']);
figure('name',figure_name)

hold on

for i = 1 : 1 : nSubresultants
   
    k = myLowerLimit + (i-1);
    
    vR1_diags = arr_R1_diag{i};
    
    vec_k = k .* ones(length(vR1_diags),1);
    
    plot(vec_k,log10(vR1_diags),'*');
    
end

vline(lowerLimit);
vline(upperLimit);

hold off

end