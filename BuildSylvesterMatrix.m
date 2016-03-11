function S = BuildSylvesterMatrix(fxy_matrix,gxy_matrix,k1,k2,alpha,theta1,theta2)
% Given two input polynomials f(x,y) and g(x,y), build the k1,k2-th
% Sylvester subresultant.
%
%       Inputs
%
% fxy_matrix : Coefficients of input polynomial f(x,y).
%
% gxy_matrix : Coefficients of input polynomial g(x,y).
%
% k1 : With respect to x.
%
% k2 : With respect to y.
%
% alpha : Optimal value of alpha.
%
% theta1 : Optimal value of theta_{1}.
%
% theta2 : Optimal value of theta_{2}.
%
%       Outputs.
%
% S : The Sylvester Subresultant S_{k_{1},k_{2}}


% Get degrees m1 and m2 of polynomial f
[r,c] = size(fxy_matrix);
m1 = r - 1;
m2 = c - 1;

% Get degrees n1 and n2 of polynomial g
[r,c] = size(gxy_matrix);
n1 = r - 1;
n2 = c - 1;

%% Preprocess polynomials f(x,y) and g(x,y)

% Preprocess polynomial f(x,y) to obtain f(w,w)
th1_mat = diag(theta1.^(0:1:m1));
th2_mat = diag(theta2.^(0:1:m2));
fww_matrix = th1_mat * fxy_matrix * th2_mat;

% Preprocess polynomial g(x,y) to obtain g(w,w)
th1_mat = diag(theta1.^(0:1:n1));
th2_mat = diag(theta2.^(0:1:n2));
gww_matrix = th1_mat * gxy_matrix * th2_mat;

% Build the partitions of the Sylvester matrix
T1 = BuildT1(fww_matrix,n1-k1,n2-k2);
T2 = BuildT1(gww_matrix,m1-k1,m2-k2);

% Build the sylvester matrix
S = [T1 alpha.*T2]; 
end