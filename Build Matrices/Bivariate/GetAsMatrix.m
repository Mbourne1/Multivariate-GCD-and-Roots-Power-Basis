function fxy_matrix = GetAsMatrix (fxy_vec, m1, m2)
%
% % Inputs
%
% fxy_vec
%
% m1 : (Int)
%
% m2 : (Int)


switch SETTINGS.VECTORISATION_METHOD
    
    case 'Version 1'
        
        fxy_matrix = GetAsMatrix_Version1(fxy_vec, m1, m2);
        
    case 'Version 2'
        
        fxy_matrix = GetAsMatrix_Version2(fxy_vec, m1, m2);
        
    otherwise
        
        error('Error')
        
end
end 
