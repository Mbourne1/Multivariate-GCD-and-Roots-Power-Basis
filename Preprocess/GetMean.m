function lambda = GetMean(fxy)
% Get the mean of the entries of f(x,y) in the matrix C(f), by a predefined
% mean method.

global SETTINGS

switch SETTINGS.MEAN_METHOD
    case 'Geometric Mean Matlab Method'
        
        lambda  = geomean(abs(nonzeros(fxy)));
    case 'None'
        lambda = 1;
        
    otherwise
        error('error')
end
