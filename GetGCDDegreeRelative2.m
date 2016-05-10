function [t1,t2] = GetGCDDegreeRelative2(fxy,gxy,m,n,t)


% Get degree structure of polynomial f(x,y)
[m1,m2] = GetDegree(fxy);

% Get degree structure of polynomial g(x,y)
[n1,n2] = GetDegree(gxy);

% %
% Produce the set of all possible t1 and t2 values
method = 'Refined';

switch method
    case 'All'
        
        % The total of t1+t2 must be between t and 2t
        k1k2_pairs = [];
        
        for t1 = min(m1,n1):-1:0;
            for t2 = min(m2,n2):-1:0;
                
                % Get new pair
                new_row = [t1 t2];
                
                % Add new pair to matrix
                k1k2_pairs = [k1k2_pairs ; new_row];
                
            end
        end
        
    case 'Refined'
        
        % Initialise the t1,t2 pair matrix
        k1k2_pairs = [];
        
        % for all possible values of t1
        for t1 = t:-1:0;
        
            % for all possible values of t2
            for t2 = t:-1:0;
            
                % Perform a series of tests
                
                % Set conditions on the degree of d(x,y)
                condM = t1 + t2 >= t;
                condA = t1 + t2 <= 2*t;
                condC = (t1 <= n1);
                condD = (t1 <= m1);
                condE = (t2 <= n2);
                condF = (t2 <= m2);
                

                % Setting bounds on polynomial v(x,y)
                condG = (n1 - t1) + (n2 - t2) <= 2*(n-t);
                condN = (n1 - t1) + (n2 - t2) >= (n-t);
                condI = n1 - t1 <= n-t;
                condJ = n2 - t2 <= n-t;
                
                % Setting bounds on polynomial u(x,y)
                condH = (m1 - t1) + (m2 - t2) <= 2*(m-t);           
                condB = (m1 - t1) + (m2 - t2) >= (m-t);
                condK = m1 - t1 <= m-t;
                condL = m2 - t2 <= m-t;
                
                
                
                
                % if all conditions are satisfied, then add the pair to the
                % list.
                if (condA && condB && condC && condD && condE && condF...
                        && condG && condH && condI && condJ && condK && condL && condM && condN)
                    
                    % Get the pair t1 t2 as a new row
                    new_row = [t1 t2];
                    % Add the new row to the pairs matrix
                    k1k2_pairs = [k1k2_pairs ; new_row];
                end
                
            end
        end
end

% %
% Remove duplicate rows of the matrix of (t1,t2) pairs
k1k2_pairs = unique(k1k2_pairs,'rows');

% Get number of pairs in the list
[nPairs,~] = size(k1k2_pairs);

display(k1k2_pairs);

% if only one [t1,t2] pair exists, then this is the degree structure.
if nPairs == 1
    fprintf('Only one possible combination of (t1,t2) \n')
    t1 = k1k2_pairs(1,1);
    t2 = k1k2_pairs(1,2);
    return
end

vMinimumSingularValues = zeros(nPairs,1);

% For each of the pairs [k1,k2]
for i = 1:1:nPairs
    
    k1 = k1k2_pairs(i,1);
    k2 = k1k2_pairs(i,2);
    
    % Build the partitions of the Sylvester matrix
    T1 = BuildT1(fxy,n1-k1,n2-k2);
    T2 = BuildT1(gxy,m1-k1,m2-k2);
    
    % Build the sylvester matrix
    Sk1k2 = [T1 T2];
    min_sing_val = min(svd(Sk1k2));
    
    vSingularValues = svd(Sk1k2);
    vMinimumSingularValues(i) = min(vSingularValues);

end

figure_name = sprintf('%s : Minimum Singular Values',mfilename)
figure('name',figure_name)
hold on
vMinimumSingularValues_log = log10(vMinimumSingularValues);
scatter3(k1k2_pairs(:,1),k1k2_pairs(:,2),vMinimumSingularValues_log)
hold off


sumk1k2 = sum(k1k2_pairs,2);
data = [k1k2_pairs sumk1k2 vMinimumSingularValues];

min_val = min(sumk1k2);
max_val = max(sumk1k2);

nValues = max_val - min_val + 1;

k1 = zeros(nValues,1);
k2 = zeros(nValues,1);
vMinimumSingularValues = zeros(nValues,1);

for i = min_val:1:max_val
    
    k = i - min_val + 1;
    
    data_filtered = data(data(:, 3) == i, :);
    
    % Get the minimum singular value of all the [k1,k2] pairs
    [val,index] = min(data_filtered(:,4));
    
    k1(k) = data_filtered(index,1);
    k2(k) = data_filtered(index,2);
    vMinimumSingularValues(k) = data_filtered(index,4);
    
end

figure_name = sprintf('%s : Singular Values',mfilename);
figure('name',figure_name)
hold on
plot(log10(vMinimumSingularValues),'-s')
hold off


% if only one value (t1+t2) 
if (nValues ==1)
   t1 = k1;
   t2 = k2;
   return
end

% Get maximum change in Singular values
[maxChange,index] = max(diff(log10(vMinimumSingularValues)));

if (maxChange < 2)
    fprintf('Insignificant Change')
    fprintf('All subresultants are either full rank or rank deficient')
    fprintf('All subresultants are full rank')
    t1 = k1(end);
    t2 = k2(end);
else
    fprintf('Significant Change')
    t1 = k1(index)
    t2 = k2(index)
end


fprintf('----------------------------------------------------------------\n')
fprintf('\n')
fprintf('The Calculated Degree of the GCD is given by \n')
fprintf('Degree of GCD wrt x : t1 = %i\n',t1)
fprintf('Degree of GCD wrt y : t2 = %i\n',t2)
fprintf('\n')
fprintf('----------------------------------------------------------------\n')

end
