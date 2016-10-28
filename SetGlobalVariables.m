function [] =  SetGlobalVariables(problem_type,ex_num,emin,emax,mean_method,bool_alpha_theta, low_rank_approx_method)

global SETTINGS


SETTINGS.PROBLEM_TYPE = problem_type;
SETTINGS.EX_NUM = ex_num;

SETTINGS.EMIN = emin;
SETTINGS.EMAX = emax;

SETTINGS.PLOT_GRAPHS = 'n';

SETTINGS.SEED = 1024;

%--------------------------------------------------------------------------
%
%       SETTINGS : PREPROCESSING
%
%


SETTINGS.MEAN_METHOD = mean_method;

SETTINGS.BOOL_ALPHA_THETA    = bool_alpha_theta;

SETTINGS.LOW_RANK_APPROXIMATION_METHOD = low_rank_approx_method;

SETTINGS.BOOL_DEGREE_METHOD = '2';



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

SETTINGS.MAX_ITERATIONS_SNTLN = 100;

SETTINGS.MAX_ERROR_SNTLN = 1e-14;


%--------------------------------------------------------------------------
%
%       SETTINGS : 
%
%
%


%
% Total
% Relative
% Both
%
SETTINGS.CALC_METHOD = 'Both';


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

