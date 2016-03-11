function dxy_matrix = GetGCDCoefficients_total(fxy_matrix,gxy_matrix,...
    uxy_matrix, vxy_matrix,...
    opt_alpha,opt_theta_1,opt_theta_2,m,n,t)
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


%% Preprocess u and v

% Preprocess u(x,y) to obtain u(w,w)
th1_mat = diag(opt_theta_1.^(0:1:m-t));
th2_mat = diag(opt_theta_2.^(0:1:m-t));

uww_matrix = th1_mat * uxy_matrix * th2_mat;

% Preprocess v(x,y) to obtain v(w,w)
th1_mat = diag(opt_theta_1.^(0:1:n-t));
th2_mat = diag(opt_theta_2.^(0:1:n-t));

vww_matrix = th1_mat * vxy_matrix * th2_mat;

%% Build Matrix C

% Build the Cauchy matrix of coefficients of u(w,w)
C1 = BuildC1_total(uww_matrix,m,t);

% Build the Cauchy matrix of coefficients of v(w,w)
C2 = BuildC1_total(vww_matrix,n,t);

% Build the RHS vector of coefficients of f and g
C = [C1;C2];

%% Build Vector f(w,w)
% Get f with thetas included

% padd f so that it is in terms of its total degree

[r,c] = size(fxy_matrix);
m1 = r -1;
m2 = c -1;
padd = zeros(m+1,m+1);
padd(1:m1+1,1:m2+1) = fxy_matrix;
fxy_matrix = padd;

[r,c] = size(gxy_matrix);
n1 = r -1;
n2 = c -1;
padd = zeros(n+1,n+1);
padd(1:n1+1,1:n2+1) = gxy_matrix;
gxy_matrix = padd;

th1_mat = diag(opt_theta_1.^(0:1:m));
th2_mat = diag(opt_theta_2.^(0:1:m));
fww_matrix = th1_mat * fxy_matrix * th2_mat;

% Get fww_matrix as a vector
fww_vec = getAsVector(fww_matrix);

% Remove the zeros associated with the polynomial by total degree
num_elements_f = nchoosek(m+2,2);
fww_vec = fww_vec(1:num_elements_f);

%% Build Vector g(w,w)

th1_mat = diag(opt_theta_1.^(0:1:n));
th2_mat = diag(opt_theta_2.^(0:1:n));
gww_matrix = th1_mat * gxy_matrix * th2_mat;

% get gww_matrix as a vector
gww_vec = getAsVector(gww_matrix);

% Remove the zeros associated with the polynomial by total degree
num_elements_g = nchoosek(n+2,2);
gww_vec = gww_vec(1:num_elements_g);

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

% Padd d(w,w) with zeros so that it can be put back into matrix form
dww_vec = [dww_vec ; zeros(nchoosek(t+2-1,2),1)];
% Arrange dw into a matrix form based on its dimensions.
dww_calc_mtrx = getAsMatrix(dww_vec,t,t);

% % Obtain d(x,y) from d(w,w)
th1_mat = diag(1./(opt_theta_1.^(0:1:t)));
th2_mat = diag(1./(opt_theta_2.^(0:1:t)));
dxy_matrix = th1_mat * dww_calc_mtrx * th2_mat;


end