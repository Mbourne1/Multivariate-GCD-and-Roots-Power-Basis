function dxy_matrix = GetGCDCoefficients_Relative_3Polys(fxy, gxy, hxy, uxy, vxy, wxy)
% Given the matrices of coefficients of f(x,y) and g(x,y), the quotient
% polynomials u(x,y) and v(x,y), and optimal values for alpha, theta_{1}
% and theta_{2}, calculate the coefficients of the GCD d(x,y), by forming
% the matrix-vector product C(u)*d = f and C(v)*d = g, and solve for d. 
%
% % Inputs
%
% [fxy, gxy, hxy] : Matrix of coefficients of polynomial f(x,y), g(x,y) and
% h(x,y)
%
% [uxy, vxy, wxy] : Matrix of coefficients of polynomial u(x,y), v(x,y) and
% w(x,y)
%
% % Outputs
%
% dxy : Coefficients of polynomial d(x,y)



% Calculate the GCD of two bivariate polynomials f(x,y) and g(x,y)

% Get degrees of polynomial f(x,y)
[m1, m2] = GetDegree(fxy);

% Get degrees of polynomial u(x,y)
[m1_t1, m2_t2] = GetDegree(uxy);

% Get degrees of polynomial d(x,y)
t1 = m1-(m1_t1);
t2 = m2-(m2_t2);

% %
% %
% Build Matrix C

% Build the Cauchy matrix of coefficients of u(\omega_{1},\omega_{2})
C1 = BuildT1_Relative(uxy,t1,t2);

% Build the Cauchy matrix of coefficients of v(\omega_{1},\omega_{2})
C2 = BuildT1_Relative(vxy,t1,t2);

% Build the Cauchy matrix of coefficients of w(\omega_{1},\omega_{2})
C3 = BuildT1_Relative(wxy,t1,t2);

% Build the RHS vector of coefficients of f and g
C = [C1; C2; C3];


% % 
% Build Vector f(w,w)

% Get fww_matrix as a vector
fxy_vec = GetAsVector(fxy);

% %
% Get gxy_matrix as a vector
gxy_vec = GetAsVector(gxy);

% Get hxy matrix as a vector
hxy_vec = GetAsVector(hxy);

% % Build the RHS vector
rhs_vec = [fxy_vec;
           gxy_vec;
           hxy_vec];

% % Calculate the solution vector

% Calculate the x vector by pinv       
x = SolveAx_b(C,rhs_vec);
dxy_vec = x;

% Arrange dw into a matrix form based on its dimensions.
dxy_matrix = GetAsMatrix(dxy_vec,t1,t2);


end