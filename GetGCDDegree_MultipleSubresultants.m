function [t] = GetGCDDegree_MultipleSubresultants(vMinimumSingularValues,degree_limits)


% Get the name of the function which called this function
[St,~] = dbstack();
calling_function = St(2).name;

% Set Global Variables.
global SETTINGS

% Get upper and lower bound on degree of GCD
lower_lim = degree_limits(1);
upper_lim = degree_limits(2);

vMinimumSingularValues = Sanitize(vMinimumSingularValues);

% Get the differences between minimum singular value S_{k} and S_{k+1}
vDeltaMinSingVal = abs(diff(log10(vMinimumSingularValues)));

% Get the maximum change (on log scale) in minimum singular values.
[max_change,index] = max((vDeltaMinSingVal));

fprintf([mfilename ' : ' calling_function ' : ' sprintf('Largest Change in Singular Values : %4.5e \n' ,abs(max_change))]);
fprintf([mfilename ' : ' calling_function ' : ' sprintf('Threshold : %4.5e \n', SETTINGS.THRESHOLD)]);

% check if the maximum change is significant

if abs(max_change) < SETTINGS.THRESHOLD
    
    
    % Get the average of the minimum singular values
    avg = mean(log10(vMinimumSingularValues));
    fprintf([mfilename ' : ' calling_function ' : ' sprintf('Average Singular Value : %4.5e \n', avg)]);
    fprintf([mfilename ' : ' calling_function ' : ' sprintf('Threshold : %4.5e \n', SETTINGS.THRESHOLD_RANK)]);
    
    if avg < SETTINGS.THRESHOLD_RANK
        % All singular values are sufficiently small for the subresultants
        % to all be rank deficient
        t = upper_lim;
        fprintf([mfilename ' : ' calling_function ' : ' sprintf('All Rank Deficient \n')]);
    else
        % All singular values are large, and indicate all subresultants are
        % full rank
        t = 0;
        fprintf([mfilename ' : ' calling_function ' : ' sprintf('All Full Rank \n')]);
    end
    
    fprintf([mfilename ' : ' calling_function ' : ' sprintf('Degree GCD : %i \n',t)]);
    
    
else
    % max_change is sufficiently large to indicate degree of GCD.
    t = lower_lim + index - 1;

    fprintf([mfilename ' : ' calling_function ' : ' 'min < Deg(GCD) < max \n'])
    fprintf([mfilename ' : ' calling_function ' : ' sprintf('Degree GCD : %i \n',t)])
    
end