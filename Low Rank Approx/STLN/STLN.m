function [fww_lr,a_gww_lr] = STLN(fww, alpha.*gww,m,n,t,t1,t2,opt_col)


switch SETTINGS.CALC_METHOD
    case 'Total'
        fww_lr = fww;
        a_gww_lr = alpha.*gww;
        
        fprintf('No STLN Developed for total degree \n')
    case 'Relative'
        % Perform STLN to obtain low rank approximation
        [fww_lr,a_gww_lr] = STLN_Relative(fww, alpha.*gww,t1,t2,opt_col);
    case 'Both'
        fprintf('No STLN developed for both total and relative degree \n')
        % Perform STLN to obtain low rank approximation
        [fww_lr,a_gww_lr] = STLN_Both(fww, alpha.*gww,m,n,t,t1,t2,opt_col);
        
end

end