function [fxy_lr,gxy_lr] = STLN(fxy,gxy,m,n,t,t1,t2,opt_col)

global SETTINGS

switch SETTINGS.CALC_METHOD
    case 'Total'
        fxy_lr = fxy;
        gxy_lr = alpha.*gxy;
        
        fprintf('No STLN Developed for total degree \n')
        
    case 'Relative'
        
        
        % Perform STLN to obtain low rank approximation
        [fxy_lr,gxy_lr] = STLN_Relative(fxy,gxy,t1,t2,opt_col);
        
    case 'Both'
        
        % Perform STLN to obtain low rank approximation
        [fxy_lr,gxy_lr] = STLN_Both(fxy,gxy,m,n,t,t1,t2,opt_col);
        
end

end