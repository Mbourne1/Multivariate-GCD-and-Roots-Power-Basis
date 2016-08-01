function [opt_col] = GetOptimalColumn(fxy,gxy,m,n,t,t1,t2)
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
% t : Total degree of polynomial d(x,y)
%
% t1 : Degree of polynomial d(x,y) with respect to x
%
% t2 : Degree of polynomial d(x,y) with respect to y

global SETTINGS

switch SETTINGS.CALC_METHOD
    case 'Total'
        opt_col = GetOptimalColumn_Total(fxy,gxy,m,n,t);
        
    case 'Relative'
        opt_col = GetOptimalColumn_Respective(fxy,gxy,t1,t2);
        
    case 'Both'
        opt_col = GetOptimalColumn_Both(fxy,gxy,m,n,t,t1,t1);
        
    otherwise
        error('err')
end


end