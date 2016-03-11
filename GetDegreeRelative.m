function [t1,t2] = GetDegreeRelative(fxy_matrix,gxy_matrix,...
    m,n,t,...
    lambda,mu,...
    opt_alpha,opt_theta_1,opt_theta_2)

global PLOT_GRAPHS

% Get degree structure of polynomial f(x,y)
[r,c] = size(fxy_matrix);
m1 = r - 1;
m2 = c - 1;

% Get degree structure of polynomial g(x,y)
[r,c] = size(gxy_matrix);
n1 = r - 1;
n2 = c - 1;

%% 
% Preprocess

% Normalise by geometric mean
fxy_matrix_working_n = fxy_matrix ./ lambda;

% Normalise by geometric mean
gxy_matrix_working_n = gxy_matrix ./ mu;

% Preprocess polynomial f(x,y) to obtain f(w,w)
th1_mat = diag(opt_theta_1.^(0:1:m1));
th2_mat = diag(opt_theta_2.^(0:1:m2));
fww_working = th1_mat * fxy_matrix_working_n * th2_mat;

% Preprocess polynomial g(x,y) to obtain g(w,w)
th1_mat = diag(opt_theta_1.^(0:1:n1));
th2_mat = diag(opt_theta_2.^(0:1:n2));
gww_working = th1_mat * gxy_matrix_working_n * th2_mat;

%%
% Produce the set of all possible t1 and t2 values
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

%%
% Remove duplicate rows of the matrix of (t1,t2) pairs
t1t2_pair_mat = unique(t1t2_pair_mat,'rows');

% If only one entry in matrix of possible t1 t2 combinations then used this
[nPairs,~] = size(t1t2_pair_mat);
if nPairs == 1
    fprintf('Only one possible combination of (t1,t2) \n')
    t1 = t1t2_pair_mat(1,1);
    t2 = t1t2_pair_mat(1,2);
    return
end
%%

% Initialise matrix
mSingularValues = [];

% For each row in the constructed matrix of (k1,k2) pairs
for i = 1:1:nPairs
    
    k1 = t1t2_pair_mat(i,1);
    k2 = t1t2_pair_mat(i,2);
    
    % Build the partitions of the Sylvester matrix
    T1 = BuildT1(fww_working,n1-k1,n2-k2);
    T2 = BuildT1(gww_working,m1-k1,m2-k2);
    
    % Build the sylvester matrix
    Sk1k2 = [T1 opt_alpha.*T2];
    min_sing_val = min(svd(Sk1k2));
    
    try
        mSingularValues = [mSingularValues ; k1 k2 log10(min_sing_val)];
        z(k1+1,k2+1) = log10(min_sing_val);
    catch
        mSingularValues = [mSingularValues ; k1 k2 0];
        z(k1+1,k2+1) = 0;
    end
    
    
end

max_k1 = max(mSingularValues(:,1));
max_k2 = max(mSingularValues(:,2));

switch method
    case 'All'
        y = 0:min(m2,n2)+1;
        x = 0:min(m1,n1)+1;
    case 'Refined'
        x = 0:max_k1 +2;
        y = 0:max_k2 +2;
end

% Add a row and column of zeros to the z matrix
[r,c] = size(z);
z = [z zeros(r,1); zeros(1,c+1)];

%%
% Plot 3d surface
switch PLOT_GRAPHS
    case 'y'
        [x,y] = meshgrid(x,y);
        figure('name','Calculating (t1,t2)')
        hold on
        s1 = surf(x,y,z');
        title('Minimum Singular Values \lambda_{i} of S_{k_{1},k_{2}}')
        xlabel('k_{1}')
        ylabel('k_{2}')
        zlabel('log_{10} Min Sing Val')
        xlim([0,max_k1+3])
        ylim([0,max_k2+3])
        alpha(s1,0.5)
        xlabel('t_{1}')
        ylabel('t_{2}')
        hold off
        
        
        %%
        % Plot 3d data points
        figure('name','3dplot')
        hold on
        title('Minimum Singular Values in S_{t_{1},t_{2}}')
        xlabel('t_{1}')
        ylabel('t_{2}')
        scatter3(mSingularValues(:,1),mSingularValues(:,2),mSingularValues(:,3),'filled')
        grid('on')
        hold off
        
    case 'n'
    otherwise
        error('plot_graphs is either y or n')
end
%%
[rows_z,cols_z] = size(z);

% Get change in x component (k1) and a row of zeros
if rows_z == 1
    delta_z_x = 0;
else
    delta_z_x = [diff(z,1); zeros(1,cols_z)];
end

% Get change in y component (k2) and a zero col

% If only one row in z
if cols_z == 1
    delta_z_y = 0;
else
    delta_z_y = [diff(z,1,2) zeros(rows_z,1)];
end


% Get change in z component
delta_zz = delta_z_x + delta_z_y;

% Get all the changes in z excluding the zeros
delta_zz(delta_zz~=0);

criterion = max((abs((delta_zz(:)))));

global THRESHOLD

if log10(criterion) < THRESHOLD
    fprintf('Value below threshold \n')
    fprintf('All subresultants are either full rank or rank deficient')
    fprintf('Degree of GCD is either 0 or min(m,n)')
    t1 = min(m1,n1);
    t2 = min(m2,n2);
    disp(t1)
    disp(t2)
    
    return
end

% Get the [i,j] entry which has maximum change to [i+1,j] and [i,j+1]
[~, idx] = max(delta_zz(:));
[x, y] = ind2sub(size(delta_zz),idx);
% Set the values t1 and t2
t1 = x-1;
t2 = y-1;


% % If only one row is returned, then take t1 and t2
[r,~] = size(mSingularValues);
if r == 1
    t1 = mSingularValues(1,1);
    t2 = mSingularValues(1,2);
    return;
end

end
