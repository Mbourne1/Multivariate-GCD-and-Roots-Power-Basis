function Sk = BuildSylvesterMatrix(fxy,gxy,m,n,k,k1,k2)



global SETTINGS


switch SETTINGS.DEGREE_METHOD
    
    case 'Total' 
        Sk = BuildSylvesterMatrix_Total(fxy,gxy,m,n,k);
    case 'Relative'
        Sk = BuildSylvesterMatrix_Relative(fxy,gxy,k1,k2);
    case 'Both'
        Sk = BuildSylvesterMatrix_Both(fxy,gxy,m,n,k,k1,k2);
        
end

end