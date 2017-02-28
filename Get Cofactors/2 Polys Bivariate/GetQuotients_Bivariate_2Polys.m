function [uxy,vxy] =  GetQuotients_Bivariate_2Polys(fxy, gxy, m, n, t, t1, t2)
%
%
% % Inputs
%
% [fxy, gxy] : Coefficients of polynomials f(x,y) and g(x,y)
%
% [m, n] : Total degree of f(x,y) and g(x,y)
%
% t : Total degree of common divisor
%
% [t1, t2] : Relative degree of common divisor
% 
% % Outputs
%
% [uxy, vxy] : Coefficients of polynomials u(x,y) and v(x,y)


global SETTINGS

% % Get coefficients of the quotients u(x,y) and v(x,y)
switch SETTINGS.DEGREE_METHOD
    case 'Total'
        [uxy,vxy] = ...
            GetQuotients_Total_Bivariate_2Polys(fxy, gxy, m, n, t);
        
    case 'Relative'
        [uxy,vxy] = ...
            GetQuotients_Relative_Bivariate_2Polys(fxy, gxy, t1, t2);
        
    case 'Both'
        [uxy,vxy] = ...
            GetQuotients_Both_Bivariate_2Polys(fxy, gxy, m, n, t, t1, t2);
        
    otherwise
        error([mfilename ': Calculation method must be total, relative or both'])
end