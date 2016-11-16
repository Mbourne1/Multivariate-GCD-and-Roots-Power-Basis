function [fxy_o, gxy_o, dxy_o, uxy_o, vxy_o, t, t1, t2] = o_gcd_mymethod(fxy, gxy, m,n, limits_t)
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
% fxy : Matrix of coefficients of polynomial f(x,y)
%
% gxy : Matrix of coefficients of polynomial g(x,y)
%
% m : Total degree of polynomial f(x,y)
%
% n : Total degree of polynomial g(x,y)
%
% % Outputs
%
% uxy : Calculated coefficients of u(x,y)
%
% vxy : Calculated coefficients of v(x,y)
%
% dxy : Calculated coefficients of d(x,y)

% % Initialise the global variables

% Initialise global settings
global SETTINGS

% % Preprocessing

[lambda,mu,alpha,th1,th2] = Preprocess_Relative(fxy,gxy);


% Normalise f(x,y) by geometric means
fxy_n = fxy./lambda;
gxy_n = gxy./mu;

% Get f(w,w) from f(x,y)
fww = GetWithThetas(fxy_n,th1,th2);

% Get g(w,w) from g(x,y)
gww = GetWithThetas(gxy_n,th1,th2);

% %
% %
% %
% Get the total degree of the GCD d(x,y)

% Get degree using total degree method
LineBreakMedium();

t_new = GetGCDDegree_Total_WithLimits(fww,alpha.*gww,m,n,limits_t);

LineBreakMedium();


t = t_new;

% If t = 0, then polynomials are coprime, and end function.
if t == 0
    uxy = fxy;
    vxy = gxy;
    dxy = 1;
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
        [t] = GetGCDDegree_Total_WithLimits(fww,alpha.*gww,m,n,limits_t);
        t1 = 1000;
        t2 = 1000;
        %[t1,t2] = GetGCDDegree_Relative(fww,alpha*gww,m,n,t);
        
    case 'Relative'
        
        t = 1000;
        [t1,t2] = GetGCDDegree_Relative(fww,alpha.*gww);
        
    case 'Both'
        [t] = GetGCDDegree_Total_WithLimits(fww,alpha.*gww,m,n,limits_t);
        
        [t1,t2] = GetGCDDegree_Relative_Given_t(fww,alpha*gww,m,n,t);
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
idx_col = GetOptimalColumn(fww,alpha.*gww,m,n,t,t1,t2);


% % 
% %
% %
% Perform low rank approximation method
[fxy_lr, gxy_lr, uxy_lr, vxy_lr, alpha_lr, th1_lr, th2_lr] = LowRankApprox(fxy_n,gxy_n,alpha,th1,th2,m,n,t,t1,t2,idx_col);



[fxy_lra,gxy_lra,uxy_lra,vxy_lra,dxy_lra,alpha_lra,th1_lra,th2_lra] = APF(fxy_lr, gxy_lr, uxy_lr, vxy_lr, m, n, t, alpha_lr, th1_lr, th2_lr);



% % 
% Remove thetas from u(w,w), v(w,w) and d(w,w) to obtain u(x,y) v(x,y) and
% d(x,y)
fxy_o = fxy_lra;
gxy_o = gxy_lra;
dxy_o = dxy_lra;
uxy_o = uxy_lra;
vxy_o = vxy_lra;


end
