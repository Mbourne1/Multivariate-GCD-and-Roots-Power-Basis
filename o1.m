function [uxy, vxy, dxy, t, t1, t2] = o1(fxy, gxy, m, n)
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
% Inputs.
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
% Outputs
%
% uxy : Calculated coefficients of u(x,y)
%
% vxy : Calculated coefficients of v(x,y)
%
% dxy : Calculated coefficients of d(x,y)

% % Initialise the global variables


global BOOL_PREPROC
global PLOT_GRAPHS
global LOW_RANK_APPROXIMATION_METHOD

% % Preprocessing
[lambda,mu,alpha,theta1,theta2] = Preprocess(fxy,gxy);

% Normalise f(x,y) by geometric means
fxy_n = fxy./lambda;
gxy_n = gxy./mu;

% Get f(w,w) from f(x,y)
fww = GetWithThetas(fxy_n,theta1,theta2);

% Get g(w,w) from g(x,y)
gww = GetWithThetas(gxy_n,theta1,theta2);

% Build the 0-th Subresultant of the preprocessed polynomials.
S_Preproc = BuildSylvesterMatrix(fww,alpha.*gww,0,0);

% % Get the Degree of the GCD d(x,y)

% Get degree using total degree method
[t] = GetGCDDegreeTotal(fxy,gxy,m,n,...
    lambda,mu,...
    alpha, theta1, theta2);

if t == 0
    uxy = fxy;
    vxy = gxy;
    dxy = 1;
    t = 0;
    t1 = 0;
    t2 = 0;
    return
end

fprintf('The calculated total degree is : %i',t)
fprintf('\n')

% Given the total degree, obtain t_{1} and t_{2}
[t1,t2] = GetGCDDegreeRelative(fxy,gxy,m,n,t,...
    lambda,mu,...
    alpha,theta1,theta2);


% % Get optimal column for removal from S_{t_{1},t_{2}}
opt_col = GetOptimalColumn(fww,alpha.*gww,t1,t2);

% Perform low rank approximation method

% Note. Use the notation 'lr' for low rank approximation.

switch LOW_RANK_APPROXIMATION_METHOD
    case 'Standard STLN'
        
        % Get preprocessed polynomials
        fww = GetWithThetas(fxy_n,theta1,theta2);
        gww = GetWithThetas(gxy_n,theta1,theta2);
        
        % Perform STLN to obtain low rank approximation
        [fww_lr,a_gww_lr] = STLN(fww, alpha.*gww,m,n,t1,t2,opt_col);
        
        % Remove alpha from alpha.*g(w,w)
        gww_lr = a_gww_lr./ alpha;
        
        % Get f(x,y) from f(w,w)
        fxy_lr = GetWithoutThetas(fww_lr,theta1,theta2);
        
        % Get g(x,y) from g(w,w)
        gxy_lr = GetWithoutThetas(gww_lr,theta1,theta2);
        
        fxy_n = fxy_lr;
        gxy_n = gxy_lr;
        
        S_LowRankApprox = BuildSylvesterMatrix(fww_lr,a_gww_lr,0,0);
        
    case 'Standard SNTLN'
        
        % Get the SNTLN of the Sylvester matrix
        [fxy_lr,gxy_lr,alpha_lr,theta1_lr,theta2_lr,~] = ...
            SNTLN(fxy_n,gxy_n,alpha,theta1,theta2,t1,t2,opt_col);
        
        
        % Update f(x,y) g(x,y) theta1 theta2 and alpha to their new values post
        % SNTLN
        fxy_n = fxy_lr;
        gxy_n = gxy_lr;
        
        theta1 = theta1_lr;
        theta2 = theta2_lr;
        alpha = alpha_lr;
        
        % Get f(w,w) from f(x,y)
        fww_lr = GetWithThetas(fxy_n,theta1,theta2);
        gww_lr = GetWithThetas(gxy_n,theta1,theta2);
        
        S_LowRankApprox = BuildSylvesterMatrix(fww_lr,alpha_lr.*gww_lr,0,0);
        
        
        
    case 'None'
        
    otherwise
        error('LOW_RANK_APPROXIMATION_METHOD must be set to either (Standard STLN) or (Standard SNTLN)')
end

% % Plot some graphs

% Build Subresultant of noisy unprocessed polynomials
S_Unproc = BuildSylvesterMatrix(fxy,gxy,0,0);
vSingularValues_S_Unproc = svd(S_Unproc);
vSingularValues_S_Unproc  = normalise(vSingularValues_S_Unproc);

switch PLOT_GRAPHS
    case 'y'
        % Plot singular values of fxy_working and fxy_output
        figure('name','Singular Values of S with and without SNTLN')
        hold on
        
        switch LOW_RANK_APPROXIMATION_METHOD
            case {'Standard SNTLN', 'Standard STLN'}
                vSingularValues_S_LowRank = svd(S_LowRankApprox);
                vSingularValues_S_LowRank = normalise(vSingularValues_S_LowRank);
                label = 'S(f(x,y),g(x,y)) SNTLN';
                plot(log10(vSingularValues_S_LowRank),'-s','DisplayName',label)
            case {'None'}
            otherwise
                error('err')
                
        end
        
        switch BOOL_PREPROC
            case 'y'
                vSingularValues_S_Preproc = svd(S_Preproc);
                vSingularValues_S_Preproc = normalise(vSingularValues_S_Preproc);
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

% % Get quotients u(x,y) and v(x,y)
calc_method = 'relative';

switch calc_method
    case 'total'
        [uxy,vxy] = ...
            GetQuotients_total(fxy_n,gxy_n,m,n,t,alpha,theta1,theta2);
    case 'relative'
        [uxy,vxy] = ...
            GetQuotients(fxy_n,gxy_n,t1,t2,alpha,theta1,theta2);
    otherwise
        error('calc method is either total or relative')
end


% % Get the GCD d(x,y)
switch calc_method
    case 'total'
        dxy = ...
            GetGCDCoefficients_total(fxy_n,gxy_n,uxy,vxy,alpha, theta1, theta2,m,n,t);
    case 'relative'
        dxy = ...
            GetGCDCoefficients(fxy_n,gxy_n,uxy,vxy,alpha, theta1, theta2);
    otherwise
        error('calc method is either total or relative')
end



end
