function arr_hxy = Deconvolve_Bivariate_Batch(arr_fxy, vDeg_t_fxy, vDeg_x_fxy, vDeg_y_fxy)
% Performs a bivariate deconvolution of polynomials f(x,y) and g(x,y) to
% obtain h(x,y)
%
% % Inputs
%
% arr_fxy : Array of polynomials f_{i}(x,y)
%
% vDeg_t_fxy : Total degree of polynomials f_{i}(x,y)
%
% vDeg_x_fxy : Degree of polynomials f_{i}(x,y) with respect to x
%
% vDeg_y_fxy : Degree of polynomials f_{i}(x,y) with respect to y


global SETTINGS
switch SETTINGS.DEGREE_METHOD
    case 'Total'
        
        arr_hxy = Deconvolve_Bivariate_Batch_Total(arr_fxy,vDeg_t_fxy);
                
    case 'Relative'
        
        arr_hxy = Deconvolve_Bivariate_Batch_Respective(arr_fxy);
        
    case 'Both'
        
        arr_hxy = Deconvolve_Bivariate_Batch_Both(arr_fxy, vDeg_t_fxy, vDeg_x_fxy, vDeg_y_fxy);
        
    otherwise
        
        error('err')
end

end