function [] = o_roots_mymethod(fxy_matrix,m)
%%

%                   Begin root finding algorithm

% Set the iteration number
ite_num = 1;

% Set the first entry of q to be the input polynomial f(x,y)
qx{1} = fxy_matrix;
vDegree_qx(1) = m;

% Get dimensions of polynomial f(x,y)
[m1,m2] = GetDegree(qx{ite_num});

% Set the iteration condition to true
iteration_condition = true;

% Whilst the most recently calculated GCD has a degree greater than
% zero. ie is not a constant, perform a gcd calculation on it and its
% derivative.
while iteration_condition
    
    fprintf('GCD calculation wrt x iteration : %i \n\n',ite_num)
    
    % Set inputs for gcd calculation
    f = qx{ite_num};

    % Get the derivative of f(x,y) with respect to x.
    g = Differentiate_wrt_x(qx{ite_num});
    
    % Get the total degree of f(x,y)
    m = vDegree_qx(ite_num);
    
    % Get the total degree of g(x,y)
    n = m - 1;
       
    % GCD is only a scalar with respect to x so set equal to g(x,y).
    [~,~,dxy_calc,t,t1,t2] = o1(f,g,m,n);
    
    % increment the iteration number
    ite_num = ite_num + 1;
    
    % Assign the GCD to be the newest member of the q array.
    qx{ite_num} = dxy_calc;
    
    % Store the degrees
    vDegree1_qx(ite_num) = t1;
    vDegree2_qx(ite_num) = t2;
    vDegree_qx(ite_num) = t;
    
   
    if t1 == 0
        iteration_condition = false;
    else
        iteration_condition = true;
    end
end

%% Obtain the series h_{i}
% Each h_{x}(i) is obtained by the deconvolution of q_{x}(i) and q_{x}(i+1)

% Get number of polynomials in the series of polynomials q_{i}
[~,num_entries_qx] = size(qx);

% For each pair of q_{x}(i) and q_{x}(i+1)
for i = 1 : 1 : num_entries_qx - 1
    
    % Perform the deconvolution
    hx{i} = Deconvolve_bivariate(qx{i},qx{i+1});
    display(hx{i})
    
    % Get the degree structure of h_{x}(i)
    deg1_hx(i) = vDegree1_qx(i) - vDegree1_qx(i+1);
    deg2_hx(i) = vDegree2_qx(i) - vDegree2_qx(i+1);
    
end



%% Obtain the series of polynomials w_{x}{i}
% Each w_{x}(i) is obtained by the deconvolution of h_{x}(i) and h_{x}(i+1)

% Get number of polynomials in the array of h_{x}
[~,num_entries_hx] = size(hx);

if num_entries_hx > 1
    
    % For each pair of h_{x}(i) and h_{x}(i+1) perform deconvolution
    for i = 1:1:num_entries_hx-1
        
        % Peform deconvolutions to obtain w_{x}
        
        wx{i}    = Deconvolve_bivariate(hx{i},hx{i+1});
        display(wx{i})
        
        % Get the degree structure of each w_{x}(i)
        deg1_wx(i) = deg1_hx(i) - deg1_hx(i+1);
        deg2_wx(i) = deg2_hx(i) - deg2_hx(i+1);
        
    end
    
    % Set the final w_{x}(i+1) to be equal to h_{x}(i+1)
    wx{i+1} = hx{i+1};
    
    % Set the final degree structure
    deg1_wx(i+1) = deg1_hx(i+1)
    deg2_wx(i+1) = deg2_hx(i+1)
    
else
    
    % Number of h_{i} is equal to one, so no deconvolutions to be performed
    wx{1} = hx{1};
    display(wx{1})
    % Get the degree structure of h_{x,i}
    deg1_wx(1) = deg1_hx(1);
    deg2_wx(1) = deg2_hx(1);
    
end


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now calculate w_{i} in terms of y

fprintf('################################################################\n')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Stage 2 - Perform GCD calculations for derivatives with respect to y


% Set the iteration number
ite_num = 1;

% Set qy to be equal to f(x,y)
qy{ite_num} = fxy_matrix;
deg_total_qy(1) = vDegree_qx(1);

% Get dimensions of polynomial f(x,y)
deg1_qy(1) = m1;
deg2_qy(1) = m2;

% Set the iteration condition to be true
iteration_condition = true;

% Whilst the most recently calculated GCD has a degree greater than
% zero. ie is not a constant, perform a gcd calculation on it and its
% derivative.
while iteration_condition
    
    fprintf('GCD calculation wrt y iteration : %i \n\n',ite_num)
    
    % Get inputs for GCD calculation
    f = qy{ite_num};
    
    % Get the derivative of f(x,y) with respect to y.
    g = Differentiate_wrt_y(qy{ite_num});
    
    % Get the total degree of f(x,y)
    m = deg_total_qy(ite_num);
    
    % Get the total degree of g(x,y)
    n = m - 1;
    
    % Get the GCD of f(x,y) and g(x,y)
    [~,~,dxy_calc,t,t1,t2] = o1(f,g,m,n);
    
    % increment the iteration number
    ite_num = ite_num + 1;
    
    % Assign the GCD to be the newest member of the q array.
    qy{ite_num} = dxy_calc;
    
    % Get degree structure of q_{y}(i)
    deg_total_qy(ite_num) = t;
    deg1_qy(ite_num) = t1;
    deg2_qy(ite_num) = t2;
    
    
    % Check the iteration condition.
    if t2 <= 0 
        iteration_condition = false;
    else
        iteration_condition = true;
    end
end

%% Obtain the series h_{i}

% get number of elements in the series of polynomials q_{i}
[~,c] = size(qy);

for i = 1:1:c-1
    % Perform deconvolution
    hy{i} = Deconvolve_bivariate(qy{i},qy{i+1});
    
    % Get the degree structure of h_{x}(i)
    deg1_hy(i) = deg1_qy(i) - deg1_qy(i+1);
    deg2_hy(i) = deg2_qy(i) - deg2_qy(i+1);
end

%% obtain the series w_{i} for the
[~,c] = size(hy);

if c > 1
    for i = 1:1:c-1
        wy{i}    =  Deconvolve_bivariate(hy{i},hy{i+1});
        
        deg1_wy(i) = deg1_hy(i) - deg1_hy(i+1);
        deg2_wy(i) = deg2_hy(i) - deg2_hy(i+1);
    end
    wy{i+1} = hy{i+1};
    deg1_wy(i+1) = deg1_hy(i+1);
    deg2_wy(i+1) = deg2_hy(i+1);
else
    wy{1} = hy{1};
    deg1_wy(1) = deg1_hy(1);
    deg2_wy(1) = deg2_hy(1);
end



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   Perform a series of GCD calculations on the w_{x,i}s

% get number of w_{x,i}
[~,c] = size(wx);
for i = 1:1:c
    [~,cols] = size(wx{i});
    if cols >1
        
        fxy_matrix_n = wx{i};
        gxy_matrix_n = wy{i};
        opt_alpha = 1;
        opt_theta_1 = 1;
        opt_theta_2 = 1;
        
        [~,c1] = size(fxy_matrix_n);
        [r2,~] = size(gxy_matrix_n);
        
        % The degree of polynomial f(x,y) with respect to y, gives the
        % degree t2.
        t2 = c1 - 1;
        
        % The degree of polynomial g(x,y) with respcet to x, gives the
        % degree t1
        t1 = r2-1;
        
        % Get quotients
        [uxy_matrix_clc,vxy_matrix_clc] = GetQuotients(fxy_matrix_n,gxy_matrix_n,t1,t2,opt_alpha,opt_theta_1,opt_theta_2);
        
        % Get the GCD dxy
        dxy_matrix_clc = GetGCDCoefficients(fxy_matrix_n,gxy_matrix_n,uxy_matrix_clc,vxy_matrix_clc,opt_alpha, opt_theta_1, opt_theta_2);
        
        % Overwrite wx and wy with new values
        wxy{i} = dxy_matrix_clc;
        wx{i} = uxy_matrix_clc;
        wy_new = Deconvolve_bivariate(wy{i},dxy_matrix_clc);
        wy{i} = vxy_matrix_clc;
    end
end

display(wx{1})

end