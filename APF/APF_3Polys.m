function [fxy_lr, gxy_lr, hxy_lr, uxy_lr, vxy_lr, wxy_lr, dxy_lr, alpha_lr, th1_lr, th2_lr ] ...
    = APF_3Polys(fxy, gxy, hxy, uxy, vxy, wxy, m, n, o, t, alpha, th1, th2)
%
% % Inputs
%
% [fxy, gxy, hxy]
%
% [uxy, vxy, wxy]
%
% [m, n, o]
%
% alpha :
%
% th1 :
%
% th2 :
%
% % Outputs
%
% [fxy_lr, gxy_lr, hxy_lr] :
%
% [uxy_lr, vxy_lr, wxy_lr] :
%
% dxy_lr :
%
% alpha_lr :
%
% [th1_lr, th2_lr] :


global SETTINGS



switch SETTINGS.APF_METHOD
    case 'None'
        
        % Get f(\omega_{1},\omega_{2}) from f(x,y) 
        fww = GetWithThetas(fxy,th1,th2);
        % Get g(\omega_{1},\omega_{2}) from g(x,y)
        gww = GetWithThetas(gxy,th1,th2);
        % Get h(\omega_{1},\omega_{2}) from h(x,y)
        hww = GetWithThetas(hxy,th1,th2);
        
        % Get u(\omega_{1},\omega_{2}) from u(x,y)
        uww = GetWithThetas(uxy,th1,th2);
        % Get v(\omega_{1},\omega_{2}) from v(x,y)
        vww = GetWithThetas(vxy,th1,th2);
        % Get w(\omega_{1},\omega_{2}) from w(x,y)
        www = GetWithThetas(wxy,th1,th2);
        
        % %
        % %
        % Get the coefficients of the GCD d(x,y)
        [dww_lr] = GetGCDCoefficients_3Polys(fww, alpha.*gww, hww, uww, vww, www, m, n, o, t);
        
        dxy_lr = GetWithoutThetas(dww_lr,th1,th2);
        
        
        fxy_lr = fxy;
        gxy_lr = gxy;
        hxy_lr = hxy;
        
        uxy_lr = uxy;
        vxy_lr = vxy;
        wxy_lr = wxy;
        
        alpha_lr = alpha;
        th1_lr = th1;
        th2_lr = th2;
        
    case 'Standard Nonlinear APF'
        
        error([mfilename ' : Method not developed']);
        
    case 'Standard Linear APF'
        
        error([mfilename ' : Method not developed']);
        
    otherwise 
        error([mfilename ' : Error'])
        
end


end