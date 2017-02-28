function [fxy_o, gxy_o, hxy_o, dxy_o, uxy_o, vxy_o, wxy_o, t, t1, t2] = ...
    o_gcd_mymethod_3Polys(fxy, gxy, hxy, m, n, o, limits_t)
% o_gcd_mymethod(fxy, gxy, m, n, limits_t)

% Given two bivariate polynomials, return the GCD d(x,y) and the coprime
% polynomials u(x,y) and v(x,y) where
%
%
% f/d = u
% g/d = v
%
% fv-gu = 0
%
% The polynomials are given matrix form as follows
%
%       1   y   y^2     ...
%        ___________
%   1   |___|___|___|   ...
%   x   |___|___|___|   ...
%   x^2 |___|___|___|   ...
%   ...    .  .  .
%
% % Inputs.
%
%
% [fxy, gxy, hxy] : Matrices of coefficients of polynomials f(x,y), g(x,y)
% and h(x,y)
%
% [m, n, o] : Total degree of polynomial f(x,y) g(x,y) and h(x,y)
%
%
% % Outputs
%
% [uxy, vxy, wxy] : Calculated coefficients of u(x,y) v(x,y) and w(x,y)
%
% dxy : Calculated coefficients of d(x,y)

% % Initialise the global variables

% Initialise global settings
global SETTINGS

% % Preprocessing

%[GM_fx, GM_gx, GM_hx, alpha, th1, th2] = Preprocess_Relative_3Polys(fxy,gxy,hxy);
% % Temporary
GM_fx = 1;
GM_gx = 1;
GM_hx = 1;
alpha = 1;
th1 = 1;
th2 = 1;


% Normalise f(x,y) by geometric means
fxy_n = fxy./GM_fx;
gxy_n = gxy./GM_gx;
hxy_n = hxy./GM_hx;

% Get f(\omega_{1},\omega_{2}) from f(x,y)
fww = GetWithThetas(fxy_n,th1,th2);

% Get g(\omega_{1},\omega_{2}) from g(x,y)
gww = GetWithThetas(gxy_n,th1,th2);

% Get h(\omega_{1},\omega_{2}) from h(x,y)
hww = GetWithThetas(hxy_n,th1,th2);


% %
% %
% %
% Get the total degree of the GCD d(x,y)

% Get degree using total degree method
LineBreakMedium();

t_new = GetGCDDegree_Total_WithLimits_Bivariate_3Polys(fww, alpha.*gww, hww, m, n, o, limits_t);

LineBreakMedium();


t = t_new;

% If t = 0, then polynomials are coprime, and end function.
if t == 0
    fxy_o = fxy;
    gxy_o = gxy;
    hxy_o = hxy;
    
    uxy_o = fxy;
    vxy_o = gxy;
    wxy_o = hxy;
    
    dxy_o = 1;
    
    t = 0;
    t1 = 0;
    t2 = 0;
    
    return
end

% Print out the total degree of the GCD
fprintf([mfilename ' : ' sprintf('The calculated total degree is : %i \n',t)]);
LineBreakMedium();




% %
% %
% %
% Given the total degree, compute relative degrees t_{1} and t_{2}

switch SETTINGS.DEGREE_METHOD
    case 'Total'
        % Must compute t1 and t2 so do so with relative method
        [t] = GetGCDDegree_Total_WithLimits_Bivariate_3Polys(fww, alpha.*gww, hww, ...
            m, n, o, limits_t);
        t1 = 1000;
        t2 = 1000;
        
        %[t1,t2] = GetGCDDegree_Relative(fww,alpha*gww,m,n,t);
        
    case 'Relative'
        
        t = 1000;
        [t1,t2] = GetGCDDegree_Relative_Bivariate_3Polys(fww, alpha.*gww, hww);
        
    case 'Both'
        
        [t] = GetGCDDegree_Total_WithLimits_Bivariate_3Polys(fww, alpha.*gww, hww,...
            m, n, o, limits_t);
        
        [t1,t2] = GetGCDDegree_Relative_Given_t_Bivariate_3Polys(fww, alpha*gww, hww,...
            m, n, o, t);
        
    otherwise
        error('err')
end



LineBreakMedium();
fprintf([mfilename ' : ' sprintf('The calculated relative degree is : t_{1} = %i, t_{2} = %i \n',t1,t2)])
LineBreakMedium();



% %
% %
% %
% Get optimal column for removal from S_{t_{1},t_{2}}
idx_col = GetOptimalColumn_3Polys(fww, alpha.*gww, hww, m, n, o, t, t1, t2);


% % 
% %
% %
% Perform low rank approximation method
[fxy_lr, gxy_lr, hxy_lr, uxy_lr, vxy_lr, wxy_lr, alpha_lr, th1_lr, th2_lr] = ...
    LowRankApprox_Bivariate_3Polys(fxy_n, gxy_n, hxy_n, alpha, th1, th2, m, n, o, t, t1, t2, idx_col);



[fxy_lra, gxy_lra, hxy_lra, uxy_lra, vxy_lra, wxy_lra, dxy_lra, alpha_lra, th1_lra, th2_lra] = ...
    APF_Bivariate_3Polys(fxy_lr, gxy_lr, hxy_lr, uxy_lr, vxy_lr, wxy_lr, m, n, o, t, alpha_lr, th1_lr, th2_lr);



% % 
% Remove thetas from u(w,w), v(w,w) and d(w,w) to obtain u(x,y) v(x,y) and
% d(x,y)
fxy_o = fxy_lra;
gxy_o = gxy_lra;
hxy_o = hxy_lra;

dxy_o = dxy_lra;

uxy_o = uxy_lra;
vxy_o = vxy_lra;
wxy_o = wxy_lra;

end
