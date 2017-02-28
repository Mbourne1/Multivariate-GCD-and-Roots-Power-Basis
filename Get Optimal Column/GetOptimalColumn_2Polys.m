function [idx_col] = GetOptimalColumn_2Polys(fxy, gxy, m, n, k, k1, k2)
% Get the optimal column c_{k} of the Sylvester subresultant matrix S(f,g)
% which when removed from S(f,g) gives minimal residual in the equation
% A_{k}x = c_{k} when solving for x.
%
% % Inputs
%
% fxy : Coefficients of polynomial f(x,y)
%
% gxy : Coefficients of polynomial g(x,y)
%
% m : Total degree of polynomial f(x,y)
% 
% n : Total degree of polynomial g(x,y)
%
% k : Total degree of polynomial d(x,y)
%
% k1 : Degree of polynomial d(x,y) with respect to x
%
% k2 : Degree of polynomial d(x,y) with respect to y
%
% % Outputs
%
% idx_col : Index of column to be removed from the Sylvester subresultant
% matrix.

global SETTINGS

switch SETTINGS.DEGREE_METHOD
    case 'Total'
        
        % Build the kth Sylvester matrix
        Sk = BuildT_Total_Bivariate_2Polys(fxy, gxy, m, n, k);
        
        % Get optimal column for removal
        idx_col = GetOptimalColumn_Total(Sk);
        
    case 'Relative'
        
        % Build the kth Sylvester matrix
        Sk1k2 = BuildT_Relative_Bivariate_2Polys(fxy, gxy, k1, k2);
        
        % Get optimal column for removal
        idx_col = GetOptimalColumn_Relative(Sk1k2);
        
    case 'Both'
        
        % Build the kth Sylvester matrix
        Skk1k2 = BuildT_Both_Bivariate_2Polys(fxy, gxy, m, n, k, k1, k2);
        
        % Get Optimal column for removal
        idx_col = GetOptimalColumn_Both(Skk1k2);
        
    otherwise
        error('err')
end


end