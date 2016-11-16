function [] =  SetGlobalVariables(problem_type, ex_num, emin, emax, ...
    mean_method, bool_alpha_theta, low_rank_approx_method, apf_method, degree_method)

global SETTINGS

SETTINGS.PLOT_GRAPHS = 'y';

% Problem Type : 'GCD' or 'Roots'
SETTINGS.PROBLEM_TYPE = problem_type;


%-------------------------------------------------------------------------
%
% Example Settings
%
%

% Example Number : String
SETTINGS.EX_NUM = ex_num;

% Noise upper and lower limit
SETTINGS.EMIN = emin;
SETTINGS.EMAX = emax;

SETTINGS.SEED = 1024;

%--------------------------------------------------------------------------


% Degree_method
%       Total
%       Relative
%       Both

SETTINGS.DEGREE_METHOD = degree_method;

%--------------------------------------------------------------------------
%
%       SETTINGS : PREPROCESSING
%
%


SETTINGS.MEAN_METHOD = mean_method;

SETTINGS.BOOL_ALPHA_THETA    = bool_alpha_theta;


%-------------------------------------------------------------------------
%
%       SETTINGS : DEGREE COMPUTATION
%
%
SETTINGS.THRESHOLD = 2;
SETTINGS.THRESHOLD_RANK = 2;

%--------------------------------------------------------------------------
%
%       SETTINGS : LOW RANK APPROXIMATION
%
%

% Low Rank Approximation Method
%       'Standard STLN'
%       'Standard SNTLN'
%       'None'
SETTINGS.LOW_RANK_APPROXIMATION_METHOD = low_rank_approx_method;

SETTINGS.MAX_ITERATIONS_SNTLN = 50;

SETTINGS.MAX_ERROR_SNTLN = 1e-10;

%--------------------------------------------------------------------------
%
%       SETTINGS : APPROXIMATE FACTORISATION SETTINGS
%
%
SETTINGS.APF_METHOD = apf_method;

%--------------------------------------------------------------------------
%
%       SETTINGS : DECONVOLUTION
%
%


%
% Total - Use total degree
% Respective - Use degree with respect to x and y
% Both - Use total and respective combined
%
SETTINGS.DECONVOLUTION_STYLE = 'Both';

% 
% Batch - Perform batch of deconvolutions together
% Separate - Each deconvolution is separate
%
SETTINGS.DECONVOLUTION_METHOD = 'Batch';

%
% 'From Deconvolutions'
% 'From uxy'
%
SETTINGS.HXY_METHOD = 'From uxy';


end

