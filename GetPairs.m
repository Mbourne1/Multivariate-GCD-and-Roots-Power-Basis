function t1t2_pair_mat = GetPairs(m,m1,m2,n,n1,n2,t)
% Used in GetGCDDegree
% Given that the total degree has been computed, get the set of pairs of
% (k_{1},k_{2}) from which the relative degree is computed.
% Get the set of (k1,k2) pairs for computing S_{k_{1},k_{2}}

% Method is either 'All' or 'Refined'
% All - All possible combinations of t1 = 0,...,min(m1,n1) and t2 =
% 0,...,min(m2,n2)
method = 'All';

switch method
    case 'All'
        
        % The total of t1+t2 must be between t and 2t
        t1t2_pair_mat = [];
        
        for t1 = min(m1,n1):-1:0;
            for t2 = min(m2,n2):-1:0;
                
                % Get new pair
                new_row = [t1 t2];
                
                % Add new pair to matrix
                t1t2_pair_mat = [t1t2_pair_mat ; new_row];
                
            end
        end
        
    case 'Refined'
        
        % Initialise the t1,t2 pair matrix
        t1t2_pair_mat = [];
        
        % for all possible values of t1
        for t1 = t:-1:0;
            
            % for all possible values of t2
            for t2 = t:-1:0;
                
                % Perform a series of tests
                condA = n1-t1 + n2 -t2 >= n-t;
                condB = m1-t1 + m2 -t2 >= m-t;
                condC = t1 <= n1 && t1 <= m1;
                condD = t2 <= n2 && t2 <= m2;
                condE = n1 - t1 + n2 - t2 <= 2*(n-t);
                condF = m1 - t1 + m2 - t2 <= 2*(m-t);
                condG = n1 - t1 <= n-t;
                condH = n2 - t2 <= n-t;
                condI = m1 - t1 <= m-t;
                condJ = m2 - t2 <= m-t;
                condK = t1+t2 >= t;
                
                
                if condA && condB && condC && condD && condE && condF...
                        && condG && condH && condI && condJ && condK
                    
                    % Get the pair t1 t2 as a new row
                    new_row = [t1 t2];
                    % Add the new row to the pairs matrix
                    t1t2_pair_mat = [t1t2_pair_mat ; new_row];
                end
                
            end
        end
end


%
% Remove duplicate rows of the matrix of (t1,t2) pairs
t1t2_pair_mat = unique(t1t2_pair_mat,'rows');


end