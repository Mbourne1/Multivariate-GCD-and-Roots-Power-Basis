function P = BuildP1_BothDegree_STLN(m,m1,m2,n,n1,n2,t,t1,t2,idx_col)
% BuildPt_sub(m1,m2,n1,n2,opt_col,t1,t2)
%
% Build the matrix P, used in STLN function. P is a matrix which is
% obtained from the decomposition of a column vector c_{t} into a matrix
% vector product P_{t} [f;g]
% c_{t} is the column of the Sylvester matrix S_{t}(f,g).
%
%   %   Inputs
%
% m : Total degree of polynomial f(x,y)
% 
% m1 : Degree of polynomial f(x,y) with respect to x
%
% m2 : Degree of polynomial f(x,y) with respect to y
%
% n : Total degree of polynomal g(x,y)
%
% n1 : Degree of polynomial g(x,y) with respect to x
%
% n2 : Degree of polynomial g(x,y) with respect to y
%
% t : Total degree of polynomial d(x,y)
%
% t1 : Degree of GCD d(x,y) with respect to x
%
% t2 : Degree of GCD d(x,y) with respect to y
%
% idx_col : Index of optiml column for removal from S(f,g)

% Get number of nonzeros of f(x,y).
nNonZeros_fxy = GetNumNonZeros(m1,m2,m);

% Build the coefficient matrix of thetas
% mat = ones(m1+1,m2+1);

% Get the number of coefficients in a polynomial of total degree m
nNonZeros_fxy_total = nchoosek(m+2,2);
nZeros_fxy_total = nchoosek(m+1,2);
vec = [...
    ones(nNonZeros_fxy_total,1);
    zeros(nZeros_fxy_total,1)
    ];
mat = GetAsMatrix(vec,m,m);
mat = mat(1:m1+1,1:m2+1);

% pad with zeros
%num_mult_wrt_x = n1-t1;
%num_mult_wrt_y = n2-t2;

% Produce a zero matrix to fill the space
padd_mat = zeros(m1+n1-t1+1, m2+n2-t2+1);

% % from the index of the optimal column, Get the number of multiplications
% with respec to x and number with respect to y.
[i,j] = GivenCol_GetIndex_Relative(n1-t1,n2-t2,idx_col);

%
ihat = i+1;
jhat = j+1;

%
nRows_f = m1+1;
nCols_f = m2+1;

% insert the theta matrix into the zero matrix
padd_mat(ihat:i+nRows_f, jhat:j+nCols_f) = mat;

% Get the padded matrix as a vector
vec_padd_mat = GetAsVector(padd_mat);

% Diagonalise the vector.
diag_mat_vec_padd_mat = diag(vec_padd_mat);

% Remove the zero columns
diag_mat_vec_padd_mat( :, ~any(diag_mat_vec_padd_mat,1) ) = [];  %columns


P = diag_mat_vec_padd_mat;

% Remove the columns corresponding to zeros in f
P = P(:,1:nNonZeros_fxy);

% Remove the rows corresponding to zeros of f*v
nNonZeros_fv = GetNumNonZeros(m1+n1-t1,m2+n2-t2,m+n-t);

P = P(1:nNonZeros_fv,:);


end