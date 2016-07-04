function hxy = Deconvolve_Bivariate(fxy,gxy,m,n)
% Performs a bivariate deconvolution of polynomials f(x,y) and g(x,y) to
% obtain h(x,y)

global SETTINGS
switch SETTINGS.DECONVOLUTION_METHOD
    case 'Total'
        error('err')
    case 'Respective'
        hxy = Deconvolve_Bivariate_respective(fxy,gxy);
    case 'Both'
        hxy = Deconvolve_Bivariate_both(fxy,gxy,m,n);
    otherwise
        error('err')
end

end