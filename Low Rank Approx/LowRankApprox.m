function [fxy_lr,gxy_lr,alpha, th1,th2] = LowRankApprox(fxy,gxy,alpha,th1,th2,m,n,t,t1,t2,opt_col)
% Compute low rank approximation of the sylvester matrix formed from
% coefficients of f(x,y) and g(x,y). Return the modified forms of f(x,y)
% and g(x,y).
%
% % Inputs
% 
% fxy_n : Coefficients of polynomial f(x,y)
%
% gxy_n : Coefficients of polynomial g(x,y)
%
% alpha : alpha
%
% th1 : theta_{1}
%
% th2 : theta_{2}
%
% m : Total degree of f(x,y)
%
% n : Total degree of g(x,y)
% 
% t : Total degree of d(x,y)
% 
% t1 : Degree of d(x,y) with respect to x
%
% t2 : Degree of d(x,y) with respect to y
%
% opt_col : Index of optimal column for removal from S_{}(f,g)
%
% % Outputs.
%
% fxy : Coefficients of polynomial f(x,y) with added perturbations
% 
% gxy : Coefficients of polynomial g(x,y) with added perturbations
%
% alpha : Refined alpha
% 
% th1 : Refined theta_{1}
% 
% th2 : Refined theta_{2}



% Global Settings
global SETTINGS

% Note. Use the notation 'lr' for low rank approximation.

switch SETTINGS.LOW_RANK_APPROXIMATION_METHOD
    case 'Standard STLN'
        
        % Get preprocessed polynomials
        fww = GetWithThetas(fxy,th1,th2);
        gww = GetWithThetas(gxy,th1,th2);
        
        % Get Low rank approximation by STLN
        [fww_lr, a_gww_lr ] = STLN(fww,alpha.*gww,m,n,t,t1,t2,opt_col);
        
        % Remove alpha from alpha.*g(w,w)
        gww_lr = a_gww_lr./ alpha;
        
        % Get f(x,y) from f(w,w)
        fxy_lr = GetWithoutThetas(fww_lr,th1,th2);
        
        % Get g(x,y) from g(w,w)
        gxy_lr = GetWithoutThetas(gww_lr,th1,th2);
        
        
    case 'Standard SNTLN'
        
        % Get low rank approximation by SNTLN
        [fxy_lr,gxy_lr,alpha_lr,theta1_lr,theta2_lr,~] = ...
                    SNTLN(fxy,gxy,alpha,th1,th2,t1,t2,opt_col);
        
        
        
        % Update f(x,y) g(x,y) theta1 theta2 and alpha to their new values post
        % SNTLN
        fxy = fxy_lr;
        gxy = gxy_lr;
        
        th1 = theta1_lr;
        th2 = theta2_lr;
        alpha = alpha_lr;
        
        % Get f(w,w) from f(x,y)
        fww_lr = GetWithThetas(fxy,th1,th2);
        gww_lr = GetWithThetas(gxy,th1,th2);
   
    case 'None'
        fxy_lr = fxy;
        gxy_lr = gxy;
    otherwise
        error([mfilename ' : ' 'LOW_RANK_APPROXIMATION_METHOD' ...
            'must be set to either (Standard STLN) or (Standard SNTLN) or (None)'])
end