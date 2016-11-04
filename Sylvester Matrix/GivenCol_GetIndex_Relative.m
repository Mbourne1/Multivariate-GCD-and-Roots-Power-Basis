function [j1,j2] = GivenCol_GetIndex_Relative(n1,n2,col_index)
% GivenCol_GetIndex_Relative(m,n,col_index)
% Given the index of a column from a convolution matrix T_{n1,k2} 
% get the values j_{1} and j_{2} such that the coefficients of f(x,y) are
% multiplied by x^{j1}y^{j2}
%
% Inputs
%
% m :    Number of mulitplications with respect to x
%
% n :    Number of multiplications with respect to y
%
% c :    Index of column removed

% Outputs
%
% j1 : number of multiplications with respect to x
%
% j2 : number of multiplications with respect to y

% Build the matrix of i coefficients
i_matrix = ones(n1+1,n2+1);
di_mat = diag(0:1:n1);
i_matrix = di_mat * i_matrix ;
i_vec = GetAsVector(i_matrix);

% Build the matrix of j coefficients
j_matrix = ones(n1+1,n2+1);
di_mat = diag(0:1:n2);
j_matrix = j_matrix * di_mat;
j_vec = GetAsVector(j_matrix);

% Get the values j_{1} and j_{2}
j1 = i_vec(col_index);
j2 = j_vec(col_index);

end