function P1 = BuildP1_TotalDegree_STLN(m,n,idx_col,k)
% BuildPt_sub(m1,m2,n1,n2,opt_col,t1,t2)
%
% Build the matrix P, used in SNTLN function. P is a matrix which is
% obtained from the decomposition of a column vector c_{t} into a matrix
% vector product P_{t} [f;g]
% c_{t} is the column of the Sylvester matrix S_{t}(f,g).
%
% % Inputs
%
% m :    Degree of polynomial f(x,y) 
%
% n :    Degree of polynomial g(x,y) 
%
% idx_col : Optimal column for removal from S(f,g)
%
% t : Degree of GCD d(x,y) 
%
% % Outputs 
%
% P1 : Matrix P1

% Build a matrix the same as f(x,y) replacing the coefficients with ones.
mat = GetAsMatrix([...
        ones(nchoosek(m+2,2),1) ;...
        zeros(nchoosek(m+1,2),1) ...
        ],m,m);

% Produce a zero matrix to fill the space
padd_mat = zeros(m+n-k+1, m+n-k+1);
nCoefficients_fv = nchoosek(m+n-k+2,2); 

% % from the index of the optimal column, Get the number of multiplications
% with respec to x and number with respect to y.
[i,j] = GivenCol_GetIndex_Relative(n-k,n-k,idx_col);


ihat = i+1;
jhat = j+1;

nRows_f = m+1;
nCols_f = m+1;

% inser the theta matrix into the zero matrix
padd_mat(ihat:i+nRows_f, jhat:j+nCols_f) = mat;

% Get the padded matrix as a vector
vec_padd_mat = GetAsVector(padd_mat);
vec_padd_mat = vec_padd_mat(1:nCoefficients_fv);

% Diagonalise the vector.
diag_mat_vec_padd_mat = diag(vec_padd_mat);


% Remove the zero columns
diag_mat_vec_padd_mat( :, ~any(diag_mat_vec_padd_mat,1) ) = [];  %columns


P1 = diag_mat_vec_padd_mat;

end