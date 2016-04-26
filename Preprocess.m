function [lambda,mu,alpha,th1,th2] = Preprocess(fxy,gxy)
% Preprocess the polynomial f(x,y) and g(x,y), and return the geometric
% means and optimal alphas and thetas from the preprocessing.
%
% Inputs.
%
% fxy : Matrix of coefficients of polynomial f(x,y)
%
% gxy : Matrix of coefficients of polynomial g(x,y)
%
%
% Outputs.
%
%
% lambda : Geometric mean of coefficients of f(x,y)
%
% mu : Geometric mean of coefficients of g(x,y)
%
% alpha : Optimal \alpha
%
% th1 : Optimal \theta_{1}
%
% th2 : Optimal \theta_{2}

% Global variables

global BOOL_ALPHA_THETA

% Get degree of f(x,y)
[m1,m2] = GetDegree(fxy);

% Get degree of g(x,y)
[n1,n2] = GetDegree(gxy);

% Get mean of f(x,y) in entries of C_{}(f)
lambda = GetMean(fxy,n1-1,n2-1);

% Get mean of g(x,y) in entries of C_{}(g)
mu = GetMean(gxy,m1-1,m2-1);

% Normalise coefficients of f(x,y) and g(x,y) by dividing the
% coefficients by the respective means.
fxy_n = fxy ./ lambda;
gxy_n = gxy ./ mu;

switch BOOL_ALPHA_THETA
    case 'y'
        
        % Obtain optimal values of alpha, theta_{1} and theta_{2}
        [alpha, th1,th2] = OptimalAlphaAndTheta(fxy_n,gxy_n);
        
        fprintf('Optimal theta_{1}  :  %0.5e \n',th1)
        fprintf('Optimal theta_{2}  :  %0.5e \n',th2)
        fprintf('Optimal alpha :   %0.5e \n', alpha)
        
        % Get f(w,w) from f(x,y)
        fww = GetWithThetas(fxy_n,th1,th2);
        
        % Get g(w,w) from g(x,y)
        gww = GetWithThetas(gxy_n,th1,th2);
        
        % Get the coefficients of f(x,y) and f(w,w) as vectors
        v_fxy = GetAsVector(fxy);
        v_fww = GetAsVector(fww);
        
        % Get the coefficients of g(x,y) and g(w,w) as vectors
        v_gxy = GetAsVector(gxy);
        v_gww = GetAsVector(gww);
        
        % Plot the coefficients of f(x,y) and f(w,w)
        PlotCoefficients(v_fxy,v_fww,'f')
        
        % Plot the coefficients of g(x,y) and g(w,w)
        PlotCoefficients(v_gxy,v_gww,'g')
        
        
    case 'n'
        
        % Set linprog outputs to be 1
        th1 = 1;
        th2 = 1;
        alpha = 1;
        
    otherwise
        error('bool_preproc is either y or n')
end