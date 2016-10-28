function [uxy,vxy] =  GetQuotients(fxy,gxy,m,n,t,t1,t2)


global SETTINGS

% % Get coefficients of the quotients u(x,y) and v(x,y)
switch SETTINGS.CALC_METHOD
    case 'Total'
        [uxy,vxy] = ...
            GetQuotients_Total(fxy,gxy,m,n,t);
        
    case 'Relative'
        [uxy,vxy] = ...
            GetQuotients_Relative(fxy,gxy,t1,t2);
        
    case 'Both'
        [uxy,vxy] = ...
            GetQuotients_Both(fxy,gxy,m,n,t,t1,t2);
        
    otherwise
        error('calc method is either total or relative')
end