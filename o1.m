function [uxy_matrix, vxy_matrix, dxy_matrix,t,t1,t2] = o1(fxy,gxy,...
    m,n)
% Given two bivariate polynomials, return the GCD d(x,y) and the coprime
% polynomials u(x,y) and v(x,y) where
%
% o1(fxy,gxy,m,n)
%
% f/d = u
% g/d = v
%
% fv-gu = 0
%
% The polynomials are given matrix form as follows
%
%       1   y   y^2     ...
%        ___________
%   1   |___|___|___|   ...
%   x   |___|___|___|   ...
%   x^2 |___|___|___|   ...
%   ...    .  .  .
%
%   Inputs.
%
%
% fxy : Matrix of coefficients of f(x,y)
%
% gxy : Matrix of coefficients of g(x,y)
%
% m : Total degree of polynomial f(x,y)
%
% n : Degree of polynomial g(x,y)
%
% Outputs
%
% uxy_matrix : Calculated coefficients of u(x,y)
%
% vxy_matrix : Calculated coefficients of v(x,y)
%
% dxy_matrix : Calculated coefficients of d(x,y)

%% Initialise the global variables


global BOOL_PREPROC
global PLOT_GRAPHS
global LOW_RANK_APPROXIMATION_METHOD

%% Get Structure of polynomials f(x,y) and g(x,y)

% Get degree of polynomial f(x,y)
[m1,m2] = GetDegree(fxy);

% Get degree of polynomial g(x,y)
[n1,n2] = GetDegree(gxy);

%% Preprocessing
switch BOOL_PREPROC
    case 'y'
        % % Include Preprocessing
        
       
            lambda  = geomean(abs(nonzeros(fxy)));
            mu      = geomean(abs(nonzeros(gxy)));
                        
            %lambda = GetGeometricMean_Total(fxy,m)
            %mu = GetGeometricMean_Total(gxy,n)
        
        % Normalise coefficients of f(x,y) and g(x,y) by dividing the
        % coefficients by the geometric mean.
        fxy_n = fxy ./ lambda;
        gxy_n = gxy ./ mu;
        
        % Obtain optimal values of alpha, theta_{1} and theta_{2}
        [alpha, theta1,theta2] = OptimalAlphaAndTheta(fxy_n,gxy_n);
        
        
        fprintf('Optimal theta_{1}  :  %0.5e \n',theta1)
        fprintf('Optimal theta_{2}  :  %0.5e \n',theta2)
        fprintf('Optimal alpha :   %0.5e \n', alpha)
        
        % Multiply each of the coefficients by the optimal theta
        
        
        % Multiply by the coefficients a_{i,j} by theta_{1}^{i}
        % and theta_{2}^{j} to obtain f(w,w)
        fww_n = GetWithThetas(fxy_n,theta1,theta2);        
        
        % Multiply by the coefficients b_{i,j} by theta_{1}^{i}
        % and theta_{2}^{j} to obtain g(w,w)
        gww_n = GetWithThetas(gxy_n,theta1,theta2);
                
        % Get the coefficients of f(x,y) and f(w,w) as vectors
        v_fxy = (GetAsVector(fxy));
        v_fww = (GetAsVector(fww_n));
                   
        PlotCoefficients(fxy,fww_n,'f')
        PlotCoefficients(gxy,gww_n,'g')
        
        % Build subresulant of normalised and preprocessed polynomials
        S_Preproc = BuildSylvesterMatrix(fxy_n,gxy_n,0,0,alpha, theta1,theta2 );

        
    case 'n'
        % Exclude pre-processing
        
        % Set geometric means of f(x) and g(x) to be 1
        lambda = 1;
        mu = 1;
        
        % Set linprog outputs to be 1
        theta1 = 1;
        theta2 = 1;
        alpha = 1;
        
        % set the normalised f(x,y) and g(x,y)
        fxy_n = fxy ./ lambda;
        gxy_n = gxy ./ mu;
        
        % set the normalised f(w,w) and g(w,w)
        fww_n = fxy ./ lambda;
        gww_n = gxy ./ mu;
    otherwise
        error('bool_preproc is either y or n')
end



%% Get the Degree of the GCD d(x,y)

% Get degree using total degree method
[t] = GetDegreeTotal(fxy,gxy,m,n,...
    lambda,mu,...
    alpha, theta1, theta2);

if t == 0
    uxy_matrix = fxy;
    vxy_matrix = gxy;
    dxy_matrix = 1;
    t = 0;
    t1 = 0;
    t2 = 0;
    return
end

fprintf('The calculated total degree is : %i',t)
fprintf('\n')

% Given the total degree, obtain t_{1} and t_{2}
[t1,t2] = GetDegreeRelative(fxy,gxy,m,n,t,...
    lambda,mu,...
    alpha,theta1,theta2);


fprintf('----------------------------------------------------------------\n')
fprintf('\n')
fprintf('The Calculated Degree of the GCD is given by \n')
fprintf('Degree of GCD wrt x : t1 = %i\n',t1)
fprintf('Degree of GCD wrt y : t2 = %i\n',t2)
fprintf('\n')
fprintf('----------------------------------------------------------------\n')



%% Get optimal column for removal from S_{t_{1},t_{2}}
opt_col = GetOptimalColumn(fww_n,gww_n,t1,t2);

%% Get the quotient polynomials uxy and vxy

switch LOW_RANK_APPROXIMATION_METHOD
    case 'Standard STLN'
        
        % Get preprocessed polynomials
        fww = GetWithThetas(fxy_n,theta1,theta2);
        gww = GetWithThetas(gxy_n,theta1,theta2);
        
        % Perform STLN
        [fww_lr,a_gww_lr] = STLN(fww, alpha.*gww,m,n,t1,t2,opt_col);
        
        
        fxy_lr = diag(1./theta1.^(0:1:m1)) * fww_lr * diag(1./theta2.^(0:1:m2));
        gxy_lr = diag(1./theta1.^(0:1:n1)) * a_gww_lr * diag(1./theta2.^(0:1:n2)) ./ alpha;
        
        
        % Get the low rank values
        alpha_lr = alpha;
        theta1_lr = theta1;
        theta2_lr = theta2;
        
        fxy_n = fxy_lr;
        gxy_n = gxy_lr;
        
        S_LowRankApprox = BuildSylvesterMatrix(fxy_lr,gxy_lr,0,0,alpha_lr, theta1_lr, theta2_lr);

    case 'Standard SNTLN'
        
        % Get the SNTLN of the Sylvester matrix
        [fxy_lr,gxy_lr,alpha_lr,theta1_lr,theta2_lr,~] = ...
            SNTLN(fxy_n,gxy_n,alpha,theta1,theta2,t1,t2,opt_col);
        
       
        % Update f(x,y) g(x,y) theta1 theta2 and alpha to their new values post
        % SNTLN
        fxy_n = fxy_lr;
        gxy_n = gxy_lr;
        
        theta1 = theta1_lr;
        theta2 = theta2_lr;
        alpha = alpha_lr;
        S_LowRankApprox = BuildSylvesterMatrix(fxy_lr,gxy_lr,0,0,alpha_lr, theta1_lr, theta2_lr);

        
        
    case 'None'
        
    otherwise
        error('bool_SNTLN must be set to either (Standard STLN) or (Standard SNTLN)')
end


PlotSingularValues()

%% Get quotients u(x,y) and v(x,y)
calc_method = 'relative';

switch calc_method
    case 'total'
        [uxy_matrix,vxy_matrix] = ...
            GetQuotients_total(fxy_n,gxy_n,m,n,t,alpha,theta1,theta2);
    case 'relative'
        [uxy_matrix,vxy_matrix] = ...
            GetQuotients(fxy_n,gxy_n,t1,t2,alpha,theta1,theta2);
    otherwise
        error('calc method is either total or relative')
end


%% Get the GCD d(x,y)
switch calc_method
    case 'total'
        dxy_matrix = ...
            GetGCDCoefficients_total(fxy_n,gxy_n,uxy_matrix,vxy_matrix,alpha, theta1, theta2,m,n,t);
    case 'relative'
        dxy_matrix = ...
            GetGCDCoefficients(fxy_n,gxy_n,uxy_matrix,vxy_matrix,alpha, theta1, theta2);
    otherwise
        error('calc method is either total or relative')
end



end
