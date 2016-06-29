function [] = o_roots(ex_num,emin,emax,mean_method,bool_alpha_theta,low_rank_approx_method)
% o_roots(ex_num,el,mean_method,bool_alpha_theta,low_rank_approx_method)
%
% Given an example number and set of parameters, obtain the roots of the
% example polynomial, where the polynomial is in the power basis
%
% %                             Inputs
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
%       'y' - Include Preprocessing
%       'n' - Exclude Preprocessing
%
% low_rank_approx_method ('y'/'n')
%       'Standard SNTLN' : Include SNTLN
%       'None'           : Exclude SNTLN
%
% % Examples
% >> o_roots('1', 1e-10, 1e-12, 'Geometric Mean Matlab Method', 'y', 'None')

% Set the problem type, used in logging to the outputs file.
problem_type = 'Roots';

SetGlobalVariables(problem_type,ex_num,emin,emax,mean_method,bool_alpha_theta,low_rank_approx_method)

% %
% %
% Get Example

% Given the example number get the polynomial coefficients of f(x,y) and
% the total degree of f.
[fxy_exact, m] = Examples_Roots(ex_num);

% Add noise to the coefficients of polynomial f(x,y)
[fxy,~] = Noise2(fxy_exact,emin,emax);


% % 
% %
% Get roots by my method
o_roots_mymethod(fxy,m);

end

