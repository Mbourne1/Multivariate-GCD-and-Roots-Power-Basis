function plotRowNorm_degreeTotal(arr_R1_RowNorm,  myLimits, limits)
%
% % Inputs
%
% arr_R1_RowNorm : 
%
% myLimits : My Limits on degree of GCD computation
%
% limits : Actual limits on degree of GCD computation

% Get my upper and lower limit
myLowerLimit = myLimits(1);
myUpperLimit = myLimits(2);

lowerLimit = limits(1);
upperLimit = limits(2);

% Get number of subresultant matrices
nSubresultants = myUpperLimit - myLowerLimit + 1;


figure_name = sprintf([mfilename ' : ' 'Row Norms']);
figure('name',figure_name);
hold on

for i = 1:1:nSubresultants
    
    vR1RowNorm = arr_R1_RowNorm{i};
    
    k = myLowerLimit + (i-1);
    
    % Produce a vector of i
    vec_k = k.* ones(length(vR1RowNorm),1);
    
    plot(vec_k, log10(vR1RowNorm),'*');
    
end

vline(lowerLimit)
vline(upperLimit)

hold off


end