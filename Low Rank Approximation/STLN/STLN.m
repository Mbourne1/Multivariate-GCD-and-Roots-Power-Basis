function [fxy_lr,gxy_lr] = STLN(fxy,gxy,m,n,k,k1,k2,idx_col)
% STLN(fxy,gxy,m,n,k,k1,k2,idx_col)
%
%
% % Inputs.
% 
% fxy : Coefficients of polynomial f(x,y)
% 
% gxy : Coefficients of polynomial g(x,y)
%
% m : Total degree of f(x,y)
%
% n : Total degree of g(x,y)
% 
% k : Total degree of d(x,y)
%
% k1 : Degree of d(x,y) with respect to x
%
% k2 : Degree of d(x,y) with respect to y


global SETTINGS

switch SETTINGS.CALC_METHOD
    case 'Total'
        
        fxy_lr = fxy;
        gxy_lr = alpha.*gxy;
        
        fprintf('No STLN Developed for total degree \n')
        
    case 'Relative'
        
        
        % Perform STLN to obtain low rank approximation
        [fxy_lr,gxy_lr] = STLN_Relative(fxy,gxy,k1,k2,idx_col);
        
        BuildSylvesterMatrix_Respective(fxy,gxy,k1,k2);
        BuildSylvesterMatrix_Respective(fxy_lr,gxy_lr,k1,k2)
        
    case 'Both'
        
        % Perform STLN to obtain low rank approximation
        [fxy_lr,gxy_lr] = STLN_Both(fxy,gxy,m,n,k,k1,k2,idx_col);
        
end

end