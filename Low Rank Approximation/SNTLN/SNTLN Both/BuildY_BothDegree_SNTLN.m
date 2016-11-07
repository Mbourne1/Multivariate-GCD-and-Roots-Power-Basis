function Y = BuildY_BothDegree_SNTLN(x,m,m1,m2,n,n1,n2,k,k1,k2,alpha,th1,th2)
% This function builds the matrix Y_{t_{1},t_{2}}, 
% Where Y(x)*z = E_{t_{1},t_{2}}(z)*x
%
% % Inputs
%
% m : Total degeree of polynomial f(x,y) 
%
% m1 : Degree of f(x,y) with respect to x
%
% m2 : Degree of f(x,y) with respect to y
% 
% n : Total degree of polynomial g(x,y)
%
% n1 : Degree of g(x,y) with respect to x
%
% n2 : Degree of g(x,y) with respect to y
%
% k : Total degree of d(x,y)
%
% k1 : Degree of d(x,y) with respect to x
%
% k2 : Degree of d(x,y) with respect to y
%
% idx_col : Index of column removed from Sylvester subresultant matrix
%
% x_ls : Least squares solution to A_{k,k1,k2}x = c_{k,k1,k2}
%
% alpha : Optimal value of \alpha
% 
% th1 : Optimal value of \theta_{1}
%
% th2 : Optimal value of \theta_{2}



% Separate the x into x1 and x2
% The first (n1-t1+1) x (n2-t2+1) coefficients are of x1
nCoeff_x1 = GetNumNonZeros(n1-k1,n2-k2,n-k);
nCoeff_x2 = GetNumNonZeros(m1-k1,m2-k2,m-k);

nCoeff_x1_mat = (n1-k1+1) * (n2-k2+1);
nCoeff_x2_mat = (m1-k1+1) * (m2-k2+1);
nZeros_x1_mat = nCoeff_x1_mat - nCoeff_x1;
nZeros_x2_mat = nCoeff_x2_mat - nCoeff_x2;


% Get the vector of coefficients of x1
x1_vec = x(1:nCoeff_x1);

% Get the vector of coefficients of x2
x2_vec = x(nCoeff_x1+1:nCoeff_x1 + nCoeff_x2);

x1_vec = [x1_vec ; zeros(nZeros_x1_mat,1)];
x2_vec = [x2_vec ; zeros(nZeros_x2_mat,1)];

% Get x1 and x2 as matrices.
x1_mat = GetAsMatrix(x1_vec,n1-k1,n2-k2);
x2_mat = GetAsMatrix(x2_vec,m1-k1,m2-k2);

% Build the convolution matrix T_{m,m1,m2}(x1)
T1 = BuildT1_Both(x1_mat,n-k,m,m1,m2);

% Build the convolution matrix T_{n,n1,n2}(x2)
T2 = BuildT1_Both(x2_mat,m-k,n,n1,n2);



% Multiply by the thetas of f and g
th1_mat = diag(th1.^(0:1:m1));
th2_mat = diag(th2.^(0:1:m2));
fww_thetas_mat = th1_mat * ones(m1+1,m2+1) * th2_mat;
th_mat1 = GetAsVector(fww_thetas_mat);
nCoeff_f = GetNumNonZeros(m1,m2,m);
th_mat1 = th_mat1(1:nCoeff_f);


%
th1_mat = diag(th1.^(0:1:n1));
th2_mat = diag(th2.^(0:1:n2));
gww_thetas_mat = th1_mat * ones(n1+1,n2+1) * th2_mat;
th_mat2 = GetAsVector(gww_thetas_mat);
nCoeff_g = GetNumNonZeros(n1,n2,n)
th_mat2 = th_mat2(1:nCoeff_g);

% multiply

vec = [th_mat1 ; th_mat2];
th_mat = diag(vec);

% multiply by the alpha of g
Y = [T1 alpha.*T2] * th_mat;
end