function dxy_matrix = GetGCDCoefficients_Relative(fxy_matrix, gxy_matrix,...
    uxy_matrix, vxy_matrix)
% Given the matrices of coefficients of f(x,y) and g(x,y), the quotient
% polynomials u(x,y) and v(x,y), and optimal values for alpha, theta_{1}
% and theta_{2}, calculate the coefficients of the GCD d(x,y), by forming
% the matrix-vector product C(u)*d = f and C(v)*d = g, and solve for d. 
%
% % Inputs
%
% fxy : Matrix of coefficients of polynomial f(x,y)
%
% gxy : Matrix of coefficients of polynomial g(x,y)
%
% uxy : Matrix of coefficients of polynomial u(x,y)
%
% vxy : Matrix of coefficients of polynomial v(x,y)
%
% % Outputs
%
% dxy : Coefficients of polynomial d(x,y)



% Calculate the GCD of two bivariate polynomials f(x,y) and g(x,y)

% Get degrees of polynomial f(x,y)
[m1,m2] = GetDegree(fxy_matrix);

% Get degrees of polynomial u(x,y)
[m1_t1,m2_t2] = GetDegree(uxy_matrix);

% Get degrees of polynomial d(x,y)
t1 = m1-(m1_t1);
t2 = m2-(m2_t2);

% %
% %
% Build Matrix C

% Build the Cauchy matrix of coefficients of u(w,w)
C1 = BuildT1(uxy_matrix,t1,t2);

% Build the Cauchy matrix of coefficients of v(w,w)
C2 = BuildT1(vxy_matrix,t1,t2);

% Build the RHS vector of coefficients of f and g
C = [C1;C2];


% % 
% Build Vector f(w,w)

% Get fww_matrix as a vector
fxy_vec = GetAsVector(fxy_matrix);

% %
% get gxy_matrix as a vector
gxy_vec = GetAsVector(gxy_matrix);

% % Build the RHS vector
rhs_vec = [fxy_vec;
           gxy_vec];

% % Calculate the solution vector

% Calculate the x vector by pinv       
x = SolveAx_b(C,rhs_vec);
dxy_vec = x;

% Arrange dw into a matrix form based on its dimensions.
dxy_matrix = GetAsMatrix(dxy_vec,t1,t2);


end