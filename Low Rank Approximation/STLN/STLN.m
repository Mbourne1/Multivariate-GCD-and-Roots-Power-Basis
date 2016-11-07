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

switch SETTINGS.DEGREE_METHOD
    case 'Total'
        
        fprintf('No STLN Developed for total degree \n')
        
        [fxy_lr,gxy_lr] = STLN_Total(fxy,gxy,m,n,k,idx_col);
              
        
        
    case 'Relative'
        
        
        % Perform STLN to obtain low rank approximation
        [fxy_lr,gxy_lr] = STLN_Relative(fxy,gxy,k1,k2,idx_col);
        
        BuildSylvesterMatrix_Relative(fxy,gxy,k1,k2);
        BuildSylvesterMatrix_Relative(fxy_lr,gxy_lr,k1,k2);
        
    case 'Both'
        
        % Perform STLN to obtain low rank approximation
        [fxy_lr,gxy_lr] = STLN_Both(fxy,gxy,m,n,k,k1,k2,idx_col);
        
        S1 = BuildSylvesterMatrix_Both(fxy,gxy,m,n,k,k1,k2);
        S2 = BuildSylvesterMatrix_Both(fxy_lr,gxy_lr,m,n,k,k1,k2);
        vec_singValues = svd(S1);
        vec_singValues_lr = svd(S2);
        
        t1 = fxy_lr - fxy;
        t2 = gxy_lr - gxy;
        norm(t1)
        norm(t2)
        
        figure()
        plot(log10(vec_singValues))
        hold on
        plot(log10(vec_singValues_lr))
        hold off
end

end