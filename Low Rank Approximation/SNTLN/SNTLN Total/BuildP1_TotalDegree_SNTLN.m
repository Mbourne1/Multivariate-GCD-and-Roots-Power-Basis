function P = BuildP1_TotalDegree_SNTLN(m,n,k,th1,th2,idx_col)
% Build the matrix P, used in SNTLN function. P is a matrix which is
% obtained from the decomposition of a column vector c_{t} into a matrix
% vector product P_{t} [f;g]
% c_{t} is the column of the Sylvester matrix S_{t}(f,g).
%
% %   Inputs
%
% m : Total degree of f(x,y)
%
% n : Total degree of g(x,y)
%
% k : Total degree of d(x,y)
%
% th1 : Optimal value of theta_{1}
%
% th2 : Optimal value of theta_{2}
%
% idx_col : Optimal column for removal from S(f,g)


% Build the coefficient matrix of thetas
pre_thetas = diag(th1.^(0:1:m));
post_thetas = diag(th2.^(0:1:m));
mat = zeros(m+1,m+1);

%mat = pre_thetas * mat * post_thetas;

for i1 = 0:1:m
    for i2 = 0:1:m
        if (i1+i2 <= m)
            mat(i1+1,i2+1) = (th1.^i1)  * (th2.^i2);
        end
    end
end



% Produce a zero matrix to fill the space
padd_mat = zeros(m+n-k+1, m+n-k+1);


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
nCoeff_fv = nchoosek(m+n-k+2,2);
vec_padd_mat = vec_padd_mat(1:nCoeff_fv);


% Diagonalise the vector.
diag_mat_vec_padd_mat = diag(vec_padd_mat);


% Remove the zero columns
diag_mat_vec_padd_mat( :, ~any(diag_mat_vec_padd_mat,1) ) = [];  %columns


P = diag_mat_vec_padd_mat;


end