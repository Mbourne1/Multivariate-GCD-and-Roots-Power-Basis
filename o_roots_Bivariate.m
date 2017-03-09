function [] = o_roots_Bivariate(ex_num, emin, emax, mean_method, bool_alpha_theta, ...
    low_rank_approx_method, apf_method, degree_method)
% o_roots_bivar(ex_num,el,mean_method,bool_alpha_theta,low_rank_approx_method)
%
% Given an example number and set of parameters, obtain the roots of the
% example polynomial, where the polynomial is in the power basis
%
% Inputs
%
% ex_num - Example Number
%
% emin : Lower noise level
%
% emax : Upper noise level
%
% mean_method : 
%       'Geometric Mean Matlab Method'
%       'None'
%
% bool_alpha_theta ('y'/'n')
%       true - Include Preprocessing
%       false - Exclude Preprocessing
%
% low_rank_approx_method ('y'/'n')
%       'Standard SNTLN' : Include Nonlinear SNTLN
%       'Standard STLN' : Standard Linear STLN
%       'None'           : Exclude SNTLN
%
% apf_method
%       'Standard APF'
%       'None'
%
% degree_method
%       'Total'
%       'Relative'
%       'Both'
%
% % Examples
% >> o_roots_Bivariate('1', 1e-10, 1e-12, 'None', false, 'None', 'None', 'Both')
% >> o_roots_Bivariate('1', 1e-10, 1e-12, 'Geometric Mean Matlab Method', true , 'Standard STLN', 'Standard APF', 'Both')


restoredefaultpath
% %
% Add subfolders
addpath(...
    'Formatting',...
    'GCD Methods',...
    'Get GCD Coefficients',...
    'Get GCD Degree',...
    'Get Optimal Column',...
    'Plotting',...
    'RootFindingMethods',...
    'Sylvester Matrix');

addpath(genpath('APF'));
addpath(genpath('Build Matrices'));
addpath(genpath('Examples'));
addpath(genpath('Get Cofactors'));
addpath(genpath('Get GCD Degree'));
addpath(genpath('Preprocessing'));
addpath(genpath('Low Rank Approximation'));
addpath(genpath('Deconvolution'))


% Set the problem type, used in logging to the outputs file.
problem_type = 'Roots';


if emax < emin
   temp_min = emax;
   emax = emin;
   emin = temp_min;
   
end

SetGlobalVariables(problem_type, ex_num, emin, emax, mean_method, ...
    bool_alpha_theta, low_rank_approx_method, apf_method, degree_method)

% %
% %
% Get Example

% Given the example number get the polynomial coefficients of f(x,y) and
% the total degree 'm' of f(x,y).
[fxy_exact, m] = Examples_Roots_Bivariate(ex_num);

% Add noise to the coefficients of polynomial f(x,y)
[fxy,~] = AddVariableNoiseToPoly(fxy_exact, emin, emax);


% % 
% %
% Get roots by my method

root_finding_method = '3 Poly GCD';

switch root_finding_method 
    case '2 Poly GCD'
        
        o_roots_mymethod_Bivariate(fxy, m);        
        
    case '3 Poly GCD'

        o_roots_mymethod_newmethod_Bivariate(fxy, m);
        
    otherwise
        
        error([mfilename ' : Error \n'])
end

end