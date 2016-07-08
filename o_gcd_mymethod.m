function [fxy,gxy,dxy,uxy, vxy, t, t1, t2] = o_gcd_mymethod(fxy, gxy, m,n, limits_t)
% Given two bivariate polynomials, return the GCD d(x,y) and the coprime
% polynomials u(x,y) and v(x,y) where
%
% o1(fxy,gxy,m,n)
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
% fxy : Matrix of coefficients of f(x,y)
%
% gxy : Matrix of coefficients of g(x,y)
%
% m : Total degree of polynomial f(x,y)
%
% n : Degree of polynomial g(x,y)
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
[lambda,mu,alpha,th1,th2] = Preprocess(fxy,gxy);

% Normalise f(x,y) by geometric means
fxy_n = fxy./lambda;
gxy_n = gxy./mu;

% Get f(w,w) from f(x,y)
fww = GetWithThetas(fxy_n,th1,th2);

% Get g(w,w) from g(x,y)
gww = GetWithThetas(gxy_n,th1,th2);

% Build the 0-th Subresultant of the preprocessed polynomials.
S_Preproc = BuildSylvesterMatrix(fxy,gxy,0,0);

% %
% %
% %
% Get the total degree of the GCD d(x,y)

% Get degree using total degree method
%LineBreakMedium();
%t_old = GetGCDDegreeTotal(fww,alpha.*gww,m,n);

LineBreakMedium();
t_new = GetGCDDegreeTotal2(fww,alpha.*gww,m,n,limits_t);
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

%[t1,t2] = GetGCDDegreeRelative(fww,alpha.*gww,m,n,t);

[t1,t2] = GetGCDDegreeRelative_Given_t(fww,alpha*gww,m,n,t);
LineBreakMedium();
fprintf([mfilename ' : ' sprintf('The calculated relative degree is : t_{1} = %i, t_{2} = %i \n',t1,t2)])
LineBreakMedium();

% 
% Total
% Respective
% Both
% 
SETTINGS.CALC_METHOD = 'Relative';

% %
% %
% %
% Get optimal column for removal from S_{t_{1},t_{2}}

opt_col = GetOptimalColumn(fww,alpha.*gww,m,n,t,t1,t2);



% % 
% %
% %
% Perform low rank approximation method

[fxy_n,gxy_n,alpha,th1,th2] = LowRankApprox(fxy_n,gxy_n,alpha,th1,th2,t1,t2,opt_col);
fww = GetWithThetas(fxy_n,th1,th2);
gww = GetWithThetas(gxy_n,th1,th2);


% %
% %
% %
% Plot some graphs

% Build Subresultant of noisy unprocessed polynomials
S_Unproc = BuildSylvesterMatrix(fxy,gxy,0,0);
vSingularValues_S_Unproc = svd(S_Unproc);
vSingularValues_S_Unproc  = Normalise(vSingularValues_S_Unproc);

switch SETTINGS.PLOT_GRAPHS
    case 'y'
        % Plot singular values of fxy_working and fxy_output
        figure('name','Singular Values of S with and without SNTLN')
        hold on
        
        switch SETTINGS.LOW_RANK_APPROXIMATION_METHOD
            case {'Standard SNTLN', 'Standard STLN'}
                vSingularValues_S_LowRank = svd(S_LowRankApprox);
                vSingularValues_S_LowRank = Normalise(vSingularValues_S_LowRank);
                label = 'S(f(x,y),g(x,y)) SNTLN';
                plot(log10(vSingularValues_S_LowRank),'-s','DisplayName',label)
            case {'None'}
            otherwise
                error('err')
                
        end
        
        switch SETTINGS.BOOL_ALPHA_THETA
            case 'y'
                vSingularValues_S_Preproc = svd(S_Preproc);
                vSingularValues_S_Preproc = Normalise(vSingularValues_S_Preproc);
                label = 'S(f(\omega_{1},\omega_{2}),g(\omega_{1},\omega_{2})) Preprocessed';
                plot(log10(vSingularValues_S_Preproc),'-o','DisplayName',label)
            case 'n'
            otherwise
                error('err')
        end
        
        label = 'S(f(x,y),g(x,y)) Noisy';
        plot(log10(vSingularValues_S_Unproc),'-o','DisplayName',label)
        legend(gca,'show');
        ylabel('log_{10} \sigma_{i} / \sigma_{1}')
        xlabel('i')
    case 'n'
    otherwise
        error('err')
end

% % 
% %
% %
% Get Cofactor polynomials u(x,y) and v(x,y)
[uww_matrix,vww_matrix] = ...
            GetQuotients(fww,alpha.*gww,m,n,t,t1,t2);



% %
% %
% Get the coefficients of the GCD d(x,y)

[dww_matrix] = GetGCDCoefficients(fww,alpha.*gww,uww_matrix,vww_matrix,m,n,t);

dxy = GetWithoutThetas(dww_matrix,th1,th2);
uxy = GetWithoutThetas(uww_matrix,th1,th2);
vxy = GetWithoutThetas(vww_matrix,th1,th2);



end
