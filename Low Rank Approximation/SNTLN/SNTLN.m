function [fxy_lr, gxy_lr, uxy_lr, vxy_lr, alpha_lr,th1_lr,th2_lr] = ...
    SNTLN(fxy, gxy, alpha, th1, th2, m, n, t, t1, t2, idx_col)
%
% % Inputs
%
% fxy : Coefficients of polynomial f(x,y)
%
% gxy : Coefficients of polynomial g(x,y)
%
% alpha : Optimal value of \alpha
%
% th1 : Optimal value of \theta_{1}
%
% th2 : Optimal value of \theta_{2}
%
% m : Total degree of polynomial f(x,y)
%
% n : Total degree of polynomial g(x,y)
%
% t : Total degree of polynomial d(x,y)
%
% t1 : Degree of d(x,y) with respect to x
%
% t2 : Degree of d(x,y) with respcet to y
%
% idx_col : Index of the optimal column to be removed from the Sylvester
% subresultant matrix.
%
% % Outputs.
%
% fxy_lr : Coefficients of polynomial f(x,y) + \delta f(x,y)
%
% gxy_lr : Coefficients of polynomial g(x,y) + \delta g(x,y)
%
% uxy_lr : Coefficients of polynomial u(x,y) + \delta u(x,y)
%
% vxy_lr : Coefficients of polynomial g(x,y) + \delta g(x,y)
%
% alpha_lr : \alpha + \delta \alpha
%
% th1_lr : \theta_{1} + \delta \theta_{1}
%
% th1_lr : \theta_{2} + \delta \theta_{2}
%
% x_lr : The solution vector x in the problem Ax = b.

% Initialise Global Settings
global SETTINGS

% Get method to compute SNTLN
switch SETTINGS.DEGREE_METHOD
    case 'Total'
        %
        [fxy_lr, gxy_lr, uxy_lr, vxy_lr, alpha_lr, th1_lr, th2_lr] = ...
            SNTLN_Total(fxy, gxy, m, n, alpha, th1, th2, t, idx_col);
        
    case 'Relative'
        % Get the SNTLN of the Sylvester matrix
        [fxy_lr, gxy_lr, uxy_lr, vxy_lr, alpha_lr, th1_lr, th2_lr] = ...
            SNTLN_Relative(fxy, gxy, alpha, th1, th2, t1, t2, idx_col);
        
    case 'Both'
        
        % Get the SNTLN of the Sylvester matrix
        [fxy_lr, gxy_lr, uxy_lr, vxy_lr, alpha_lr, th1_lr, th2_lr] = ...
            SNTLN_Both(fxy, gxy, alpha, th1, th2, m, n, t, t1, t2, idx_col);
end

end