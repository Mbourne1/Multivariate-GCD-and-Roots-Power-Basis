function [dxy] = GetGCDCoefficients(fxy,gxy,uxy,vxy,m,n,k)
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
    
        dxy = GetGCDCoefficients_Total(fxy, gxy, uxy, vxy, m, n, k);
        
    case 'Relative'
        
        dxy = GetGCDCoefficients_Relative(fxy, gxy, uxy, vxy);
    
    case 'Both'
        
        dxy = GetGCDCoefficients_Both(fxy, gxy, uxy, vxy, m, n, k);
        
    otherwise
        
        error([mfilename ' : DEGREE_METHOD must be Total, Relative, or Both'])
        
end