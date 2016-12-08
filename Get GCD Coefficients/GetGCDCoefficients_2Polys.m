function [dxy] = GetGCDCoefficients_2Polys(fxy, gxy, uxy, vxy, m, n, k)
% GetGCDCoefficients(fxy,gxy,uxy,vxy,m,n,k)
%
% % Inputs
%
% [fxy, gxy] : Coefficients of the polynomial f(x,y) and g(x,y)
%
% [uxy, vxy] : Coefficients of the polynomial u(x,y) and v(x,y)
%
% [m, n] : Total degree of polynomial f(x,y) and g(x,y)
%
% t : Total degree of polynomial d(x,y)

global SETTINGS

switch SETTINGS.DEGREE_METHOD
    case 'Total'
    
        dxy = GetGCDCoefficients_Total_2Polys(fxy, gxy, uxy, vxy, m, n, k);
        
    case 'Relative'
        
        dxy = GetGCDCoefficients_Relative_2Polys(fxy, gxy, uxy, vxy);
    
    case 'Both'
        
        dxy = GetGCDCoefficients_Both_2Polys(fxy, gxy, uxy, vxy, m, n, k);
        
    otherwise
        
        error([mfilename ' : DEGREE_METHOD must be Total, Relative, or Both'])
        
end