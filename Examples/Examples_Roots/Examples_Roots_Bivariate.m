function [fxy_matrix, m] = Examples_Roots_Bivariate(ex_num)
% Given an example number, get a set of roots, Build the polynomial and
% output the matrix of coefficients.
%
% Inputs.
%
% ex_num : Example Number
%
%
% Outputs.
%
% fxy_matrix : Coefficient matrix of polynomial f(x,y)
%
% m : Total degree of polynomial f(x,y)

EXAMPLE_TYPE = 'From Coefficients';

switch EXAMPLE_TYPE
    case 'From Roots'
        [fxy_matrix, m] = Examples_Roots_FromRoots_Bivariate(ex_num);
        
    case 'From Coefficients'
        [fxy_matrix, m] = Examples_Roots_FromCoefficients_Bivariate(ex_num);
        
end


end


