function [] = o_roots_mymethod(fxy_matrix,M)
%%

%                   Begin root finding algorithm

% Set the iteration number
ite = 1;

% Set the first entry of q to be the input polynomial f(x,y)
fx{1} = fxy_matrix;

vDegt_fx(ite) = M(ite);

% Get dimensions of polynomial f(x,y)
[m1,m2] = GetDegree(fx{ite});

vDeg1_fx(ite) = m1;
    
vDeg2_fx(ite) = m2;

% Whilst the most recently calculated GCD has a degree greater than
% zero. ie is not a constant, perform a gcd calculation on it and its
% derivative.
while vDeg1_fx(ite) > 0
    
    fprintf('GCD Calculation Loop iteration = %i \n', ite );
    fprintf('Compute GCD of f_{%i} and derivative f_{%i}\n\n',ite,ite);   
    
    % Get the derivative of f(x,y) with respect to x.
    gxy = Differentiate_wrt_x(fx{ite});
    
    % Get the total degree of f(x,y)
    m =  vDegt_fx(ite);
    
    % Get the total degree of g(x,y)
    n =  m - 1;
       
    if ite > 1
        lower_lim = vDegt_fx(ite)-dx(ite-1);   
        upper_lim = m-1;
    else
        lower_lim = 1;
        upper_lim = m-1;
    end
    
    fprintf('Minimum degree of f_{%i}: %i \n', ite+1, lower_lim);
    fprintf('Maximum degree of f_{%i}: %i \n', ite+1, upper_lim);
    
    % GCD is only a scalar with respect to x so set equal to g(x,y).
    [fx{ite},~,fx{ite+1},~,~,t,t1,t2] = o_gcd_mymethod(fx{ite},gxy,m,n,[lower_lim, upper_lim]);

    % Get degree of d(x,y) with respect to x
    vDeg1_fx(ite+1) = t1;
    
    % Get degree of d(x,y) with respect to y
    vDeg2_fx(ite+1) = t2;
    
    % Get total degree of d(x,y)
    vDegt_fx(ite+1) = t;
        
    % Get number of distinct roots of f(ite)
    dx(ite) = vDegt_fx(ite) - vDegt_fx(ite+1);
    
    fprintf('The computed deg(GCD(f_{%i},f_{%i}) is : %i \n',ite,ite,vDegt_fx(ite+1))
    fprintf('Number of distinct roots in f_{%i} : %i \n',ite,dx(ite))
    fprintf('Degree of f_{%i} : %i \n',ite + 1, vDegt_fx(ite+1))
    
    % increment the iteration number
    ite = ite + 1;
    
end

% % Obtain the series h_{i}
% Each h_{x}(i) is obtained by the deconvolution of q_{x}(i) and q_{x}(i+1)

% Get number of polynomials in the series of polynomials q_{i}
[~,num_entries_qx] = size(fx);

% For each pair of q_{x}(i) and q_{x}(i+1)
for i = 1 : 1 : num_entries_qx - 1
    
    % Perform the deconvolution
    hx{i} = Deconvolve_bivariate(fx{i},fx{i+1});

    % Get the degree structure of h_{x}(i)
    deg1_hx(i) = vDeg1_fx(i) - vDeg1_fx(i+1);
    deg2_hx(i) = vDeg2_fx(i) - vDeg2_fx(i+1);
    
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now calculate w_{i} in terms of y
fprintf('################################################################\n')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

% Stage 2 - Perform GCD calculations for derivatives with respect to y


% Set the iteration number
ite = 1;

% Set qy to be equal to f(x,y)
fy{ite} = fxy_matrix;
vDegt_fy(ite) = vDegt_fx(1);

% Get dimensions of polynomial f(x,y)
vDeg1_fy(ite) = m1;
vDeg2_fy(ite) = m2;



% Whilst the most recently calculated GCD has a degree greater than
% zero. ie is not a constant, perform a gcd calculation on it and its
% derivative.
while vDeg2_fy(ite) > 0
    
    fprintf('GCD calculation wrt y iteration : %i \n\n',ite)
    
    % Get inputs for GCD calculation
    f = fy{ite};
    
    % Get the derivative of f(x,y) with respect to y.
    gxy = Differentiate_wrt_y(fy{ite});
    
    % Get the total degree of f(x,y)
    m = vDegt_fy(ite);
    
    % Get the total degree of g(x,y)
    n = m - 1;
    
    if ite > 1
        lower_lim = vDegt_fy(ite)-dy(ite-1);   
        upper_lim = m-1;
    else
        lower_lim = 1;
        upper_lim = m-1;
    end
    
    fprintf('Minimum degree of f_{%i}: %i \n', ite+1, lower_lim);
    fprintf('Maximum degree of f_{%i}: %i \n', ite+1, upper_lim);
    
    % Get the GCD of f(x,y) and g(x,y)
    [~,~,dxy_calc,~,~,t,t1,t2] = o_gcd_mymethod(f,gxy,m,n,[lower_lim,upper_lim]);
    
    % Assign the GCD to be the newest member of the q array.
    fy{ite+1} = dxy_calc;
    
    % Get degree structure of q_{y}(i)
    %
    vDegt_fy(ite+1) = t;
    
    %
    vDeg1_fy(ite+1) = t1;
    
    %
    vDeg2_fy(ite+1) = t2;
    
    % Get number of distinct roots of f(ite)
    dy(ite) = vDegt_fx(ite) - vDegt_fx(ite+1);
    
    fprintf('The computed deg(GCD(f_{%i},f_{%i}) is : %i \n',ite,ite,vDegt_fx(ite+1))
    fprintf('Number of distinct roots in f_{%i} : %i \n',ite,dx(ite))
    fprintf('Degree of f_{%i} : %i \n',ite + 1, vDegt_fx(ite+1))
    
    % increment the iteration number
    ite = ite + 1;
end

%% Obtain the series h_{i}

% get number of elements in the series of polynomials q_{i}
[~,c] = size(fy);

for i = 1:1:c-1
    % Perform deconvolution
    hy{i} = Deconvolve_bivariate(fy{i},fy{i+1});
    
    % Get the degree structure of h_{x}(i)
    deg1_hy(i) = vDeg1_fy(i) - vDeg1_fy(i+1);
    deg2_hy(i) = vDeg2_fy(i) - vDeg2_fy(i+1);
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