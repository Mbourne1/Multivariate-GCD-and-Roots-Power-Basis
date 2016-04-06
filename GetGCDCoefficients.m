function dxy_matrix = GetGCDCoefficients(fxy_matrix,gxy_matrix,...
    uxy_matrix, vxy_matrix,...
    opt_alpha,th1,th2)
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
% th1 : \theta_{2}
%
% th2 : \theta_{1}

% Calculate the GCD of two bivariate polynomials f(x,y) and g(x,y)

%% Get Degree structures

% Get degrees of polynomial f
[m1,m2] = GetDegree(fxy_matrix);

% Get degrees of polynomial u
[m1_t1,m2_t2] = GetDegree(uxy_matrix);

% Get degrees of polynomial d
t1 = m1-(m1_t1);
t2 = m2-(m2_t2);

%% Preprocess u and v

% Preprocess u(x,y) to obtain u(w,w)
uww_matrix = GetWithThetas(uxy_matrix,th1,th2);

% Preprocess v(x,y) to obtain v(w,w)
vww_matrix = GetWithThetas(vxy_matrix,th1,th2);

%% Build Matrix C

% Build the Cauchy matrix of coefficients of u(w,w)
C1 = BuildT1(uww_matrix,t1,t2);

% Build the Cauchy matrix of coefficients of v(w,w)
C2 = BuildT1(vww_matrix,t1,t2);

% Build the RHS vector of coefficients of f and g
C = [C1;C2];

%% Build Vector f(w,w)
% Get f with thetas included
fww_matrix = GetWithThetas(fxy_matrix,th1,th2);

% Get fww_matrix as a vector
fww_vec = GetAsVector(fww_matrix);

%% Build Vector g(w,w)
gww_matrix = GetWithThetas(gxy_matrix,th1,th2);

% get gww_matrix as a vector
gww_vec = GetAsVector(gww_matrix);

%% Build the RHS vector
rhs_vec = [fww_vec;
           opt_alpha .*gww_vec];

%% Calculate the solution vector

% Calculate the x vector by pinv       
x = SolveAx_b(C,rhs_vec);
dww_vec = x;

% Arrange dw into a matrix form based on its dimensions.
dww_calc_mtrx = GetAsMatrix(dww_vec,t1,t2);

% % Obtain d(x,y) from d(w,w)
dxy_matrix = GetWithoutThetas(dww_calc_mtrx,th1,th2);



end