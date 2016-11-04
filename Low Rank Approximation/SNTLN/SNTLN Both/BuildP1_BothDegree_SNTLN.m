function P = BuildP1_BothDegree_SNTLN(m,m1,m2,n,n1,n2,th1,th2,idx_col,k,k1,k2)
% Build the matrix P, used in SNTLN function. P is a matrix which is
% obtained from the decomposition of a column vector c_{t} into a matrix
% vector product P_{t} [f;g]
% c_{t} is the column of the Sylvester matrix S_{t}(f,g).
%
% %   Inputs
%
% m : Total degree of f(x,y)
%
% m1 :    Degree of polynomial f(x,y) with respect to x
%
% m2 :    Degree of polynomial f(x,y) with respect to y
%
% n : Total degree of g(x,y)
%
% n1 :    Degree of polynomial g(x,y) with respect to x
%
% n2 :    Degree of polynomial g(x,y) with respect to y
%
% th1 : Optimal value of theta_{1}
%
% th2 : Optimal value of theta_{2}
%
% idx_col : Optimal column for removal from S(f,g)
%
% k : Total degree of d(x,y)
%
% k1 : Degree of GCD d(x,y) with respect to x
%
% k2 : Degree of GCD d(x,y) with respect to y


% Build the coefficient matrix of thetas
pre_thetas = diag(th1.^(0:1:m1));
post_thetas = diag(th2.^(0:1:m2));
mat = zeros(m1+1,m2+1);
%mat = pre_thetas * mat * post_thetas;

for i1 = 0:1:m1
    for i2 = 0:1:m2
        if (i1+i2 <= m)
            mat(i1+1,i2+1) = (th1.^i1)  * (th2.^i2);
        end
    end
end



% Produce a zero matrix to fill the space
padd_mat = zeros(m1+n1-k1+1, m2+n2-k2+1);

% % from the index of the optimal column, Get the number of multiplications
% with respec to x and number with respect to y.
[i,j] = GivenCol_GetIndex_Relative(n1-k1,n2-k2,idx_col);


ihat = i+1;
jhat = j+1;

nRows_f = m1+1;
nCols_f = m2+1;

% inser the theta matrix into the zero matrix
padd_mat(ihat:i+nRows_f, jhat:j+nCols_f) = mat;

% Get the padded matrix as a vector
vec_padd_mat = GetAsVector(padd_mat);

nCoeff_fv = GetNumNonZeros(m1+n1-k1,m2+n2-k2,m+n-k);




% Diagonalise the vector.
diag_mat_vec_padd_mat = diag(vec_padd_mat);


% Remove the zero columns
diag_mat_vec_padd_mat( :, ~any(diag_mat_vec_padd_mat,1) ) = [];  %columns


P = diag_mat_vec_padd_mat;
P = P(1:nCoeff_fv,:);

end