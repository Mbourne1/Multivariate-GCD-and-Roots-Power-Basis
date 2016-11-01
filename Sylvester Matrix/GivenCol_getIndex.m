function [i,j] = GivenCol_GetIndex(m,n,col_index)
% GivenCol_GetIndex(m,n,col_index)
%
% Inputs
%
% m :    Number of mulitplications with respect to x
%
% n :    Number of multiplications with respect to y
%
% c :    Index of column removed

% Outputs
% i : number of multiplications with respect to x
% j : number of multiplications with respect to y

% build the matrix of i coefficients
i_matrix = ones(m+1,n+1);

% multiply the rows by 0,1,2,3,4
di_mat = diag(0:1:m);
i_matrix = di_mat * i_matrix ;

i_vec = GetAsVector(i_matrix);

% Build the matrix of j coefficients
j_matrix = ones(m+1,n+1);

di_mat = diag(0:1:n);
j_matrix = j_matrix * di_mat;

j_vec = GetAsVector(j_matrix);

i = i_vec(col_index);
j = j_vec(col_index);

end