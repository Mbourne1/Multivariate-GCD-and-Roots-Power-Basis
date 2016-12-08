function Sk = BuildSylvesterMatrix_3Polys(fxy, gxy, hxy, m, n, o, k, k1, k2)
%
% % Inputs
%
% [fxy, gxy, hxy] : Coefficients of polynomials f(x,y), g(x,y) and h(x,y)
%
% [m, n, o] : Total degrees of polynomials f(x,y), g(x,y) and h(x,y)
%
% k : Total degree of d(x,y)
%
% [k1, k2] : Relative degrees of d(x,y)

global SETTINGS


switch SETTINGS.DEGREE_METHOD
    
    case 'Total' 
        
        Sk = BuildSylvesterMatrix_Total_3Polys(fxy, gxy, hxy, m, n, o, k);
        
    case 'Relative'
        
        Sk = BuildSylvesterMatrix_Relative_3Polys(fxy, gxy, hxy, k1, k2);
        
    case 'Both'
        
        Sk = BuildSylvesterMatrix_Both_3Polys(fxy, gxy, hxy, m, n, o, k, k1, k2);
        
end

end