function [dxy] = GetGCDCoefficients_3Polys(fxy, gxy, hxy, uxy, vxy, wxy, m, n, o, k)
% GetGCDCoefficients(fxy,gxy,uxy,vxy,m,n,k)
%
% % Inputs
%
% fxy : Coefficients of the polynomial f(x,y)
%
% gxy : Coefficients of the polynomial g(x,y)
%
% uxy : Coefficients of the polynomial u(x,y)
%
% vxy : Coefficients of the polynomial v(x,y)
%
% m : Total degree of polynomial f(x,y)
%
% n : Total degree of polynomial g(x,y)
%
% t : Total degree of polynomial d(x,y)

global SETTINGS

switch SETTINGS.DEGREE_METHOD
    case 'Total'
    
        dxy = GetGCDCoefficients_Total_3Polys(fxy, gxy, hxy, uxy, vxy, wxy, m, n, o, k);
        
    case 'Relative'
        
        dxy = GetGCDCoefficients_Relative_3Polys(fxy, gxy, hxy, uxy, vxy, wxy);
    
    case 'Both'
        
        dxy = GetGCDCoefficients_Both_3Polys(fxy, gxy, hxy, uxy, vxy, wxy, m, n, o, k);
        
    otherwise
        
        error([mfilename ' : DEGREE_METHOD must be Total, Relative, or Both'])
        
end