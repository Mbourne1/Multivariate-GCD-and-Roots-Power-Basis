function dxy_matrix = GetGCDCoefficients(fxy_matrix,gxy_matrix,...
    uxy_matrix, vxy_matrix,...
    opt_alpha,opt_theta_1,opt_theta_2)
% Given the matrices of coefficients of f(x,y) and g(x,y), the quotient
% polynomials u(x,y) and v(x,y), and optimal values for alpha, theta_{1}
% and theta_{2}, calculate the coefficients of the GCD d(x,y).
%
% %         Inputs
%
% fxy_matrix
%
% gxy_matrix
%
% uxy_matrix
%
% vxy_matrix
%
% opt_alpha
%
% opt_theta_1
%
% opt_theta_2

% Calculate the GCD of two bivariate polynomials f(x,y) and g(x,y)

%% Get Degree structures

% Get degrees of polynomial f
[r,c] = size(fxy_matrix);
m1 = r - 1;
m2 = c - 1;

% Get degrees of polynomial g
[r,c] = size(gxy_matrix);
n1 = r - 1;
n2 = c - 1;

% Get degrees of polynomial u
[r,c] = size(uxy_matrix);
m1_t1 = r - 1;
m2_t2 = c - 1;

% Get degrees of polynomial v
[r,c] = size(vxy_matrix);
n1_t1 = r - 1;
n2_t2 = c - 1;

% Get degrees of polynomial d
t1 = m1-(m1_t1);
t2 = m2-(m2_t2);

%% Preprocess u and v

% Preprocess u(x,y) to obtain u(w,w)
th1_mat = diag(opt_theta_1.^(0:1:m1-t1));
th2_mat = diag(opt_theta_2.^(0:1:m2-t2));

uww_matrix = th1_mat * uxy_matrix * th2_mat;

% Preprocess v(x,y) to obtain v(w,w)
th1_mat = diag(opt_theta_1.^(0:1:n1-t1));
th2_mat = diag(opt_theta_2.^(0:1:n2-t2));

vww_matrix = th1_mat * vxy_matrix * th2_mat;

%% Build Matrix C

% Build the Cauchy matrix of coefficients of u(w,w)
C1 = BuildC1(uww_matrix,t1,t2,m1,m2);

% Build the Cauchy matrix of coefficients of v(w,w)
C2 = BuildC1(vww_matrix,t1,t2,n1,n2);

% Build the RHS vector of coefficients of f and g
C = [C1;C2];

%% Build Vector f(w,w)
% Get f with thetas included

th1_mat = diag(opt_theta_1.^(0:1:m1));
th2_mat = diag(opt_theta_2.^(0:1:m2));
fww_matrix = th1_mat * fxy_matrix * th2_mat;

% Get fww_matrix as a vector
fww_vec = GetAsVector(fww_matrix);

%% Build Vector g(w,w)

th1_mat = diag(opt_theta_1.^(0:1:n1));
th2_mat = diag(opt_theta_2.^(0:1:n2));
gww_matrix = th1_mat * gxy_matrix * th2_mat;

% get gww_matrix as a vector
gww_vec = GetAsVector(gww_matrix);

%% Build the RHS vector
rhs_vec = [fww_vec;
           opt_alpha .*gww_vec];

%% Calculate the solution vector

% Calculate the x vector by pinv       
x = pinv(C) * rhs_vec;
dww_vec = x;

% % Calculate the x vector by QR decomposition
% [~,n2] = size(C);
% [Q,R] = qr(C);
% R1 = R(1:n2,:);
% cd = Q'*rhs_vec;
% c = cd(1:n2,:);
% x_ls = R1\c;
% dww_vec = x_ls;
    
% Get the residual associated with the solution x. (Small residual implies good approximation)    
residual = pinv(C)*rhs_vec - x;

% Arrange dw into a matrix form based on its dimensions.
dww_calc_mtrx = GetAsMatrix(dww_vec,t1,t2);

% % Obtain d(x,y) from d(w,w)
th1_mat = diag(1./(opt_theta_1.^(0:1:t1)));
th2_mat = diag(1./(opt_theta_2.^(0:1:t2)));
dxy_matrix = th1_mat * dww_calc_mtrx * th2_mat;


end