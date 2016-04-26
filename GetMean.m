function lambda = GetMean(fxy,n1_k1,n2_k2)
% Get the mean of the entries of f(x,y) in the matrix C(f), by a predefined
% mean method.

global MEAN_METHOD

switch MEAN_METHOD
    case 'Geometric Mean Matlab Method'
        
        lambda  = geomean(abs(nonzeros(fxy)));
    case 'None'
        lambda = 1;
        
    otherwise
        error('err')
end
