function [] =  SetGlobalVariables(bool_preproc, low_rank_approx_method)

global BOOL_PREPROC
BOOL_PREPROC    = bool_preproc;

global LOW_RANK_APPROXIMATION_METHOD
LOW_RANK_APPROXIMATION_METHOD = low_rank_approx_method;

global bool_degreeMethod
bool_degreeMethod = '2';

global bool_remove_cols
bool_remove_cols = 'n';

global PLOT_GRAPHS
PLOT_GRAPHS = 'y';

global SEED 
SEED = 1024;

global THRESHOLD
THRESHOLD = 1;

global MAX_ITERATIONS_SNTLN
MAX_ITERATIONS_SNTLN = 100;

global MAX_ERROR_SNTLN
MAX_ERROR_SNTLN = 1e-15;


end

