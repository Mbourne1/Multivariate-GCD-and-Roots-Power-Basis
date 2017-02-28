function plotSingularValues_degreeTotal(arr_SingularValues, limits)
% 
% % Inputs
%
% arr_SingularValues
%
% limits


figure()
hold on

for i = 1:1:length(arr_SingularValues)
   
    % Get vector of singular values of S_{k}
    vSingularValues = arr_SingularValues{i};
    
    v_is = i.*ones(length(vSingularValues),1);
    
    % plot singular values
    plot(v_is, log10(vSingularValues),'*')
    
end

hold off


end