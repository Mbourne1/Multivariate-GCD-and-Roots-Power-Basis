function [] =  SetGlobalVariables(mean_method,bool_alpha_theta, low_rank_approx_method)

global SETTINGS

SETTINGS.MEAN_METHOD = mean_method;

SETTINGS.BOOL_ALPHA_THETA    = bool_alpha_theta;

SETTINGS.LOW_RANK_APPROXIMATION_METHOD = low_rank_approx_method;

SETTINGS.BOOL_DEGREE_METHOD = '2';

SETTINGS.BOOL_REMOVE_COLS = 'n';

SETTINGS.PLOT_GRAPHS = 'n';

SETTINGS.SEED = 1024;

SETTINGS.THRESHOLD = 1;

SETTINGS.MAX_ITERATIONS_SNTLN = 100;

SETTINGS.MAX_ERROR_SNTLN = 1e-15;


end

