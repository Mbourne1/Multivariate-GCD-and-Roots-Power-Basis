function Sk = BuildSylvesterMatrix_2Polys(fxy, gxy, m, n, k, k1, k2)
%
% % Inputs
%
% [fxy, gxy] : Coefficients of the polynomials f(x,y) and g(x,y) 
%
% [m, n] : Degree of f(x,y) and g(x,y)
%
% k :
%
% [k1, k2] :
%
% % Outputs



global SETTINGS


switch SETTINGS.DEGREE_METHOD
    
    case 'Total' 
        
        Sk = BuildSylvesterMatrix_Total_2Polys(fxy, gxy, m, n, k);
        
    case 'Relative'
        
        Sk = BuildSylvesterMatrix_Relative_2Polys(fxy, gxy, k1, k2);
        
    case 'Both'
        
        Sk = BuildSylvesterMatrix_Both_2Polys(fxy, gxy, m, n, k, k1, k2);
        
end

end