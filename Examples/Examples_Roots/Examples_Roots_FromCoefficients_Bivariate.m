function [fxy, m] = Examples_Roots_FromCoefficients_Bivariate(ex_num)
%
% % Inputs
%
% ex_num : Example Number
%
% % Outputs
%
% fxy : Coefficients of polynomial f(x,y)

% Add path to examples files
addpath(genpath('../Examples'))

% Get the symbolic roots and their multiplicities
[f_root_sym_mult_arr] = Roots_Examples_Bivariate(ex_num);

% Get the coefficients of the bivariate polynomial f(x,y)
[fxy, m, m1, m2] = GetCoefficientsFromSymbolicRoots(f_root_sym_mult_arr);





end