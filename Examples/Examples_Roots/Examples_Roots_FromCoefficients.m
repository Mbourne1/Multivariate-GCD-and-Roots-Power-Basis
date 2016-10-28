
function [fxy_matrix, m] = Examples_Roots_FromCoefficients(ex_num)


x = sym('x');
y = sym('y');


addpath '../Examples'

% Get the symbolic roots and their multiplicities
[f_root_sym_mult_arr] = Bivariate_Roots_Examples(ex_num);

% Get the coefficients of the bivariate polynomial f(x,y)
[fxy_matrix,m,m1,m2] = GetCoefficientsFromSymbolicRoots(f_root_sym_mult_arr);





end