function [t] = GetGCDDegree_MultipleSubresultants(vMetric, degree_limits)
%
% Given a metric - provided by one of the GetGCDDegree functions, determine
% the degree of the GCD.
%
% % Inputs
%
% vMetric : vector of values from the metric used to determine whether each
% sylvester subresultant matrix is full rank or rank deficient.
%
% degree_Limits : [lowerLimit upperLimit], index of first and last
% sylvester subresultant whose rank revealing metric is computed.
% 

% Get the name of the function which called this function
[St,~] = dbstack();
calling_function = St(2).name;

% Set Global Variables.
global SETTINGS

% Get upper and lower bound on degree of GCD
lowerLimit = degree_limits(1);
upperLimit = degree_limits(2);

% Sanitize the metric.
vMetric = Sanitize(vMetric);

% Get the differences between minimum singular value S_{k} and S_{k+1}
vDeltaMetric = abs(diff(log10(vMetric)));

% Get the maximum change (on log scale) in minimum singular values.
[max_change, index] = max((vDeltaMetric));

fprintf([mfilename ' : ' calling_function ' : ' sprintf('Largest Change in Rank Revealing Metric : %4.5e \n' ,abs(max_change))]);
fprintf([mfilename ' : ' calling_function ' : ' sprintf('Threshold : %4.5e \n', SETTINGS.THRESHOLD)]);

% check if the maximum change is significant

if abs(max_change) < SETTINGS.THRESHOLD
    
    
    % Get the average of the minimum singular values
    avg = mean(log10(vMetric));
    fprintf([mfilename ' : ' calling_function ' : ' sprintf('Average Value of Metric: %4.5e \n', avg)]);
    fprintf([mfilename ' : ' calling_function ' : ' sprintf('Threshold : %4.5e \n', SETTINGS.THRESHOLD_RANK)]);
    
    if avg < SETTINGS.THRESHOLD_RANK
        % All singular values are sufficiently small for the subresultants
        % to all be rank deficient
        t = upperLimit;
        fprintf([mfilename ' : ' calling_function ' : ' sprintf('All Subresultant matrices are rank Deficient \n')]);
    else
        % All singular values are large, and indicate all subresultants are
        % full rank
        t = 0;
        fprintf([mfilename ' : ' calling_function ' : ' sprintf('All subresultant matrices are full rank \n')]);
    end
    
    fprintf([mfilename ' : ' calling_function ' : ' sprintf('Degree of the GCD : %i \n',t)]);
    
    
else
    % max_change is sufficiently large to indicate degree of GCD.
    t = lowerLimit + index - 1;

    fprintf([mfilename ' : ' calling_function ' : ' 'min < Deg(GCD) < max \n'])
    fprintf([mfilename ' : ' calling_function ' : ' sprintf('Degree of the GCD : %i \n',t)])
    
end