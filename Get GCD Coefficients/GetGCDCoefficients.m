function [dxy] = GetGCDCoefficients(fxy,gxy,uxy,vxy,m,n,t)

global SETTINGS

switch SETTINGS.DEGREE_METHOD
    case 'Total'
        dxy = ...
            GetGCDCoefficients_Total(fxy,gxy,uxy,vxy,m,n,t);
    case 'Relative'
        dxy = ...
            GetGCDCoefficients_Relative(fxy,gxy,uxy,vxy);
    case 'Both'
        dxy = ...
            GetGCDCoefficients_Both(fxy,gxy,uxy,vxy,m,n,t);
    otherwise
        error([mfilename ' : calc method is either total or relative'])
end