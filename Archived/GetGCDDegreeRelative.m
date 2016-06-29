function [t1,t2] = GetGCDDegreeRelative(fxy,gxy,...
    m,n,t)
% Get the degree structure (t_{1} and t_{2}) of the GCD d(x,y) of the two 
% polynomials f(x,y) and g(x,y)
%
%   Inputs.
%
%   fxy : Coefficient matrix of polynomial f(x,y)
%
%   gxy : Coefficient matrix of polynomial g(x,y)
%
%   m : Total degree of polynomial f(x,y)
%
%   n : Total degree of polynomial g(x,y)
%
%   t : Total degree of GCD d(x,y)

global SETTINGS

% Get degree structure of polynomial f(x,y)
[m1,m2] = GetDegree(fxy);

% Get degree structure of polynomial g(x,y)
[n1,n2] = GetDegree(gxy);

% Produce the set of all possible t1 and t2 values
method = 'All';

% Get set of all pairs [k1,k2]
k1k1Pairs = GetPairs(m,m1,m2,n,n1,n2,t,method);

% If only one entry in matrix of possible t1 t2 combinations then used this
[nPairs,~] = size(k1k1Pairs);
if nPairs == 1
    fprintf('Only one possible combination of (t1,t2) \n')
    t1 = k1k1Pairs(1,1);
    t2 = k1k1Pairs(1,2);
    return
end
%%

% Initialise matrix
m_SingularValues = zeros(nPairs,1);

% For each row in the constructed matrix of (k1,k2) pairs
for i = 1:1:nPairs
    
    % Get the ith pair of k_{1} and k_{2} values
    k1 = k1k1Pairs(i,1);
    k2 = k1k1Pairs(i,2);
    
    % Build the partitions of the Sylvester matrix
    T1 = BuildT1(fxy,n1-k1,n2-k2);
    T2 = BuildT1(gxy,m1-k1,m2-k2);
    
    % Build the sylvester matrix
    Sk1k2 = [T1 T2];
    
    % Get the singular values of S_{k_{1},k_{2}}
    vSingularValues = svd(Sk1k2);
    
    % Get the minimal singular value
    vMinimumSingularValues_all(i) = min(vSingularValues);
    
      
    
end

% Plot the minimum singular values
PlotGraphs_DegreeRelative()




if log10(criterion) < SETTINGS.THRESHOLD
    fprintf('Value below threshold \n')
    fprintf('All subresultants are either full rank or rank deficient')
    fprintf('Degree of GCD is either 0 or min(m,n)')
    t1 = min(m1,n1);
    t2 = min(m2,n2);
    disp(t1)
    disp(t2)
    
    return
end


LineBreakMedium()
fprintf([mfilename ' : ' 'The Calculated Degree of the GCD is given by \n'])
fprintf([mfilename ' : ' sprintf('Degree of GCD wrt x : t1 = %i\n',t1)])
fprintf([mfilename ' : ' sprintf('Degree of GCD wrt y : t2 = %i\n',t2)])
LineBreakMedium()
end
