function [t1,t2] = GetGCDDegree_Relative_3Polys(fxy, gxy, hxy)
% Get the degree structure (t_{1} and t_{2}) of the GCD d(x,y) of the two
% polynomials f(x,y) and g(x,y)
%
% Inputs.
%
% fxy : Coefficient matrix of polynomial f(x,y)
%
% gxy : Coefficient matrix of polynomial g(x,y)
%
% hxy : Coefficient matrix of polynomial h(x,y)
%
% Outputs
%
% t1 : degree of d(x,y) with respect to x
% 
% t2 : degree of d(x,y) with respect to y



% Get the degree structure of polynomial f(x,y)
[m1,m2] = GetDegree(fxy);

% Get the degree structure of polynomial g(x,y)
[n1,n2] = GetDegree(gxy);

% Get the degree structure of polynomial h(x,y)
[o1,o2] = GetDegree(hxy);

% Get the set of all pairs of (k1,k2) combinations
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
    
    % Build the partitions of the Sylvester matrix
    T1 = BuildT1_Relative(fxy,n1-k1,n2-k2);
    T2 = BuildT1_Relative(fxy,o1-k1,o2-k2);
    
    T3 = BuildT1_Relative(gxy,m1-k1,m2-k2);
    T4 = BuildT1_Relative(hxy,m1-k1,m2-k2);
    
    diagonal = blkdiag(T1,T2);
    column = [T3; T4];
    % Build the sylvester matrix
    Sk1k2 = [diagonal column];
    
    % Get the singular values of S_{k_{1},k_{2}}
    vSingularValues = svd(Sk1k2);
    
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

v_k1 = zeros(nValues,1);
v_k2 = zeros(nValues,1);
vMinimumSingularValues = zeros(nValues,1);

for i = min_val:1:max_val
    
    k = i - min_val + 1;
    
    data_filtered = data(data(:, 3) == i, :);
    
    % Get the minimum singular value of all the [k1,k2] pairs
    [~,index] = min(data_filtered(:,4));
    
    v_k1(k) = data_filtered(index,1);
    v_k2(k) = data_filtered(index,2);
    vMinimumSingularValues(k) = data_filtered(index,4);
    
    
end

% Plot the minimum singular values
PlotGraphs_DegreeRelative();

% If only one value (t1+t2)
if (nValues == 1)
   t1 = v_k1;
   t2 = v_k2;
   return
end




% Get maximum change in Singular values
[maxChange,index] = max(diff(log10(vMinimumSingularValues)));

global SETTINGS

if (maxChange < SETTINGS.THRESHOLD_RANK)
    fprintf([mfilename 'Insignificant Change\n'])
    fprintf([mfilename 'All subresultants are either full rank or rank deficient\n'])
    fprintf([mfilename 'All subresultants are full rank\n'])
    t1 = v_k1(end);
    t2 = v_k2(end);
else
    fprintf([mfilename ' : ' 'Significant Change in Singular Values\n'])
    t1 = v_k1(index);
    t2 = v_k2(index);
end

fprintf([mfilename ' : ' 'The Calculated Degree of the GCD is given by \n'])
fprintf([mfilename ' : ' sprintf('Degree of GCD wrt x : t1 = %i\n',t1)])
fprintf([mfilename ' : ' sprintf('Degree of GCD wrt y : t2 = %i\n',t2)])
end
