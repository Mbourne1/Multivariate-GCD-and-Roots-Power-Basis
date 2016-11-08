function Y = BuildY_TotalDegree_SNTLN(x,m,n,k,alpha,th1,th2)
% This function builds the matrix Y_{t_{1},t_{2}}, 
% Where Y(x)*z = E_{t_{1},t_{2}}(z)*x
%
% % Inputs
%
% x : Least squares solution to A_{k,k1,k2}x = c_{k,k1,k2} with zero
% inserted.
%
% m : Total degeree of polynomial f(x,y) 
% 
% n : Total degree of polynomial g(x,y)
%
% k : Total degree of d(x,y)
%
% alpha : Optimal value of \alpha
% 
% th1 : Optimal value of \theta_{1}
%
% th2 : Optimal value of \theta_{2}
%
% % Outputs.
%
% Y : Matrix Y.


% Separate the x into x1 and x2
% The first (n1-t1+1) x (n2-t2+1) coefficients are of x1
nNonZeros_x1 = nchoosek(n-k+2,2);
nNonZeros_x2 = nchoosek(m-k+2,2);

nZeros_x1_mat = nchoosek(n-k+1,2);
nZeros_x2_mat = nchoosek(m-k+1,2);

% Get the number of non-zeros in f(x,y) and g(x,y)
nNonZeros_fxy = nchoosek(m+2,2);
nNonZeros_gxy = nchoosek(n+2,2);

% Get the vector of coefficients of x1
x1_vec = x(1:nNonZeros_x1);

% Get the vector of coefficients of x2
x2_vec = x(nNonZeros_x1+1:nNonZeros_x1 + nNonZeros_x2);

x1_vec = [x1_vec ; zeros(nZeros_x1_mat,1)];
x2_vec = [x2_vec ; zeros(nZeros_x2_mat,1)];

% Get x1 and x2 as matrices.
x1_mat = GetAsMatrix(x1_vec,n-k,n-k);
x2_mat = GetAsMatrix(x2_vec,m-k,m-k);

% Build the convolution matrix T_{m,m1,m2}(x1)
T_x1 = BuildT1_Total(x1_mat,n-k,m);

% Build the convolution matrix T_{n,n1,n2}(x2)
T_x2 = BuildT1_Total(x2_mat,m-k,n);

% Multiply the two convolution matrices T(x1) T(x2) by the thetas 
% corresponding to f(w,w) and g(w,w) so [T(x1) T(x2)]* [f;g]

% Get thetas corresponding to the coefficients of polynomial f(x,y)
th1_mat = diag(th1.^(0:1:m));
th2_mat = diag(th2.^(0:1:m));
fww_thetas_mat = th1_mat * ones(m+1,m+1) * th2_mat;
th_mat1 = GetAsVector(fww_thetas_mat);
th_mat1 = th_mat1(1:nNonZeros_fxy);


% Get thetas corresponding to the coefficients of polynomial g(x,y)
th1_mat = diag(th1.^(0:1:n));
th2_mat = diag(th2.^(0:1:n));
gww_thetas_mat = th1_mat * ones(n+1,n+1) * th2_mat;
th_mat2 = GetAsVector(gww_thetas_mat);
th_mat2 = th_mat2(1:nNonZeros_gxy);

% Get thetas as a diagonal matrix
vec = [th_mat1 ; th_mat2];
th_mat = diag(vec);

% Multiply the second partiton by alpha, and multiply the matrix T(x1)
% T(x2) by the theta matrix.
Y = [T_x1 alpha.*T_x2] * th_mat;

end