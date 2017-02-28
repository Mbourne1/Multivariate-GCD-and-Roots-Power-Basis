function [fxy_lr, gxy_lr, uxy_lr, vxy_lr] = STLN(fxy,gxy,m,n,k,k1,k2,idx_col)
% STLN(fxy,gxy,m,n,k,k1,k2,idx_col)
%
%
% % Inputs.
%
% [fxy, gxy] : Coefficients of polynomial f(x,y) and g(x,y)
%
% [m, n] : Total degree of f(x,y) and g(x,y)
%
% k : Total degree of d(x,y)
%
% k1 : Degree of d(x,y) with respect to x
%
% k2 : Degree of d(x,y) with respect to y
%
% % Outputs.
%
% fxy_lr :
%
% gxy_lr :
%
% uxy_lr :
%
% vxy_lr :

global SETTINGS

switch SETTINGS.DEGREE_METHOD
    case 'Total'
        
        
        
        [fxy_lr, gxy_lr, uxy_lr, vxy_lr] = STLN_Total(fxy, gxy, m, n, k, idx_col);
        
        
    case 'Relative'
        
        
        % Perform STLN to obtain low rank approximation
        [fxy_lr, gxy_lr, uxy_lr, vxy_lr] = STLN_Relative(fxy,gxy,k1,k2,idx_col);
        
        
        BuildT_Relative_Bivariate_2Polys(fxy, gxy, k1, k2);
        BuildT_Relative_Bivariate_2Polys(fxy_lr, gxy_lr, k1, k2);
        
    case 'Both'
        
        % Perform STLN to obtain low rank approximation
        [fxy_lr, gxy_lr, uxy_lr, vxy_lr] = STLN_Both(fxy, gxy, m, n, k, k1, k2, idx_col);
        
        
end


end