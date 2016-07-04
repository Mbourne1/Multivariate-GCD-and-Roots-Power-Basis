function arr_hxy = Deconvolve_Bivariate_Batch(arr_fxy,vDeg_fxy)
% Performs a bivariate deconvolution of polynomials f(x,y) and g(x,y) to
% obtain h(x,y)

global SETTINGS
switch SETTINGS.DECONVOLUTION_STYLE
    case 'Total'
        
        error('err')
        
    case 'Respective'
        
        arr_hxy = Deconvolve_Bivariate_Batch_Respective(arr_fxy);
        
    case 'Both'
        
        arr_hxy = Deconvolve_Bivariate_Batch_Both(arr_fxy,vDeg_fxy);
        
    otherwise
        
        error('err')
end

end