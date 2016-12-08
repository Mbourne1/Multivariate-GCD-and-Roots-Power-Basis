function [t1,t2] = GetGCDDegree_Relative_Given_t_3Polys(fxy, gxy, hxy, m, n, o, t)
% Get the degree structure (t_{1} and t_{2}) of the GCD d(x,y) of the two
% polynomials f(x,y) and g(x,y)
%
% % Inputs.
%
% [fxy, gxy, hxy] : Coefficient matrices of polynomials f(x,y), g(x,y) and
% h(x,y)
%
% [m, n, o] : Total degree of polynomial f(x,y), g(x,y) and h(x,y)
%
% t : Total degree of GCD d(x,y)
%
% % Outputs
%
% t1 : degree of d(x,y) with respect to x
% 
% t2 : degree of d(x,y) with respect to y

% Note this file differs from 'GetGCDDegreeRelative' since we make use of
% the computed value of 't' given t, for any k1,k2 pair, we know the
% location of certain zeros (if they exist in the lower right triangle) for
% polynomials u(x,y) and v(x,y). We remove the corresponding columns of
% S(f,g)

% Get the degree structure of polynomial f(x,y)
[m1, m2] = GetDegree(fxy);

% Get the degree structure of polynomial g(x,y)
[n1, n2] = GetDegree(gxy);

% Get the degree structure of polynomial h(x,y0
[o1, o2] = GetDegree(hxy);

% Get the set of all pairs of (k1,k2) combinations
%k1k2Pairs = GetPairs_Refined(m,m1,m2,n,n1,n2,t);
k1k2Pairs = GetPairs_All_3Polys(m1,m2,n1,n2,o1,o2);

% Get number of pairs in the list
[nPairs,~] = size(k1k2Pairs);

% if only one (k1,k2) pair exists, then (t1,t2) = (k1,k2)
if nPairs == 1
    fprintf([mfilename ' : ' sprintf('Only one possible combination of (t1,t2) \n')])
    t1 = k1k2Pairs(1,1);
    t2 = k1k2Pairs(1,2);
    return
end

vMinimumSingularValues_all = zeros(nPairs,1);

% For each of the pairs [k1,k2]
for i = 1:1:nPairs

    % Get the ith pair of (k1,k2)
    k1 = k1k2Pairs(i,1);
    k2 = k1k2Pairs(i,2);
    
    % Build the Sylvester matrix S_{k,k1,k2}
    Skk1k2 = BuildSylvesterMatrix_Both_3Polys(fxy, gxy, hxy, m, n, o, t, k1, k2);
    
    % Get the singular values of S_{k_{1},k_{2}}
    vSingularValues = svd(Skk1k2);
    
    % Get the minimal singular value
    vMinimumSingularValues_all(i) = min(vSingularValues);

end


% Sort all k1 k2 pairs by their total k_{total} =  k1 + k2. Get the minimum of the 
% minimum singular values for each total.

% Get the sum of k1 + k2 for all (k1,k2) pairs
sumk1k2 = sum(k1k2Pairs,2);

% Create a matrix of data
data = [k1k2Pairs sumk1k2 vMinimumSingularValues_all];

% Get the minimum sum
min_val = min(sumk1k2);

% Get the maximum sum
max_val = max(sumk1k2);

% Get the number of possible sum values.
nValues = max_val - min_val + 1;

k1 = zeros(nValues,1);
k2 = zeros(nValues,1);
vMinimumSingularValues = zeros(nValues,1);

for i = min_val:1:max_val
    
    k = i - min_val + 1;
    
    data_filtered = data(data(:, 3) == i, :);
    
    % Get the minimum singular value of all the [k1,k2] pairs
    [~,index] = min(data_filtered(:,4));
    
    k1(k) = data_filtered(index,1);
    k2(k) = data_filtered(index,2);
    vMinimumSingularValues(k) = data_filtered(index,4);
    
    
end

% Plot the minimum singular values
PlotGraphs_DegreeRelative();

% If only one value (t1+t2)
if (nValues == 1)
   t1 = k1;
   t2 = k2;
   return
end




% Get maximum change in Singular values
[maxChange,index] = max(diff(log10(vMinimumSingularValues)));

[t] = GetGCDDegree_MultipleSubresultants(vMinimumSingularValues,[min_val,max_val]);

global SETTINGS

if (maxChange < SETTINGS.THRESHOLD_RANK)
    fprintf([mfilename 'Insignificant Change\n'])
    fprintf([mfilename 'All subresultants are either full rank or rank deficient\n'])
    fprintf([mfilename 'All subresultants are full rank\n'])
    t1 = k1(end);
    t2 = k2(end);
else
    fprintf([mfilename ' : ' 'Significant Change in Singular Values\n'])
    t1 = k1(index);
    t2 = k2(index);
end

fprintf([mfilename ' : ' 'The Calculated Degree of the GCD is given by \n'])
fprintf([mfilename ' : ' sprintf('Degree of GCD wrt x : t1 = %i\n',t1)])
fprintf([mfilename ' : ' sprintf('Degree of GCD wrt y : t2 = %i\n',t2)])
end
