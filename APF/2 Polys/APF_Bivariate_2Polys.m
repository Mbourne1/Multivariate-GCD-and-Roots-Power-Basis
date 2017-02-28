function [fxy_lr,gxy_lr,uxy_lr,vxy_lr,dxy_lr,alpha_lr,th1_lr,th2_lr ] ...
    = APF_Bivariate_2Polys(fxy, gxy, uxy, vxy, m, n, t, alpha, th1, th2)
%
% % Inputs
%
% [fxy, gxy] : Coefficients of polynomials f(x,y) and g(x,y)
%
% [uxy, vxy] : Coefficients of polynomials u(x,y) and v(x,y)
%
% [m, n] : Total degree of f(x,y) and g(x,y)
%
% t : Total degree of d(x,y)
%
% alpha :
%
% th1 :
%
% th2 :

global SETTINGS



switch SETTINGS.APF_METHOD
    case 'None'
        
        % Get f(w,w) from f(x,y) and g(w,w) from g(x,y)
        fww = GetWithThetas(fxy,th1,th2);
        gww = GetWithThetas(gxy,th1,th2);
        uww = GetWithThetas(uxy,th1,th2);
        vww = GetWithThetas(vxy,th1,th2);
        
        % %
        % %
        % Get the coefficients of the GCD d(x,y)
        [dww_lr] = GetGCDCoefficients_Bivariate_2Polys(fww,alpha.*gww,uww,vww,m,n,t);
        
        dxy_lr = GetWithoutThetas(dww_lr,th1,th2);
        
        fxy_lr = fxy;
        gxy_lr = gxy;
        
        uxy_lr = uxy;
        vxy_lr = vxy;
        
        alpha_lr = alpha;
        th1_lr = th1;
        th2_lr = th2;
        
    case 'Standard Nonlinear APF'
        
        error([mfilename ' : Method not developed']);
        
    case 'Standard Linear APF'
        
        error([mfilename ' : Method not developed']);
    otherwise 
        error([mfilename ' : Error APF Method is either None, Standard Nonlinear APF or Standard Linear APF '])
        
end


end