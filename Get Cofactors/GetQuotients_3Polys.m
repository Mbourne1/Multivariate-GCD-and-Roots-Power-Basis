function [uxy, vxy, wxy] =  GetQuotients_3Polys(fxy, gxy, hxy, m, n, o, t, t1, t2)
%
% % Inputs
%
% [fxy, gxy, hxy] : 
%
% [m, n, o]
%
% t
%
% [t1, t2] 


global SETTINGS

% % Get coefficients of the quotients u(x,y) and v(x,y)
switch SETTINGS.DEGREE_METHOD
    case 'Total'
        [uxy, vxy, wxy] = ...
            GetQuotients_Total_3Polys(fxy, gxy, hxy, m, n, o, t);
        
    case 'Relative'
        [uxy, vxy, wxy] = ...
            GetQuotients_Relative_3Polys(fxy, gxy, hxy, t1, t2);
        
    case 'Both'
        [uxy, vxy, wxy] = ...
            GetQuotients_Both_3Polys(fxy, gxy, hxy, m, n, o, t, t1, t2);
        
    otherwise
        error('calc method is either total or relative')
end