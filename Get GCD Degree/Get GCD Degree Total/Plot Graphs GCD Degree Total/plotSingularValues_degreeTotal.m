function plotSingularValues_degreeTotal(arr_SingularValues, myLimits_t, limits_t)
% 
% % Inputs
%
% arr_SingularValues
%
% myLimits_t :
%
% limits_t :


myLowerLimit = myLimits_t(1);
myUpperLimit = myLimits_t(2);

lowerLimit = limits_t(1);
upperLimit = limits_t(2);

nSubresultants = myUpperLimit - myLowerLimit + 1;

figure_name = sprintf([mfilename ' : ' 'Plotting Singular Values']);
figure('name', figure_name)
hold on


for i = 1 : 1 : nSubresultants
   
    k = myLowerLimit + (i-1);
    
    % Get vector of singular values of S_{k}
    vSingularValues = arr_SingularValues{i};
    
    v_ks = k.*ones(length(vSingularValues),1);
    
    % plot singular values
    plot(v_ks, log10(vSingularValues),'*')
    
end

vline(lowerLimit);
vline(upperLimit);

hold off


end