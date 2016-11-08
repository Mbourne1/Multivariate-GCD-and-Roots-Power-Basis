function dxy_matrix = GetGCDCoefficients_Both(fxy_matrix,gxy_matrix,...
    uxy_matrix, vxy_matrix,m,n,k)
% Given the matrices of coefficients of f(x,y) and g(x,y), the quotient
% polynomials u(x,y) and v(x,y), and optimal values for alpha, theta_{1}
% and theta_{2}, calculate the coefficients of the GCD d(x,y).
%
% %         Inputs
%
% fxy_matrix : Coefficients of polynomial f(x,y)
%
% gxy_matrix : Coefficients of polynomial g(x,y)
%
% uxy_matrix : Coefficients of polynomial u(x,y)
%
% vxy_matrix : Coefficients of polynomial v(x,y)
%
% m : Total degree of f(x,y)
% 
% n : Total degree of g(x,y)
% 
% k : Total degree of d(x,y)

% Calculate the GCD of two bivariate polynomials f(x,y) and g(x,y)

% %
% %
% Get Degree structures

% Get degrees of polynomial f(x,y)
[m1,m2] = GetDegree(fxy_matrix);

% Get degree of polynomial g(x,y)
[n1,n2] = GetDegree(gxy_matrix);

% Get degrees of polynomial u(x,y)
[m1_k1,m2_k2] = GetDegree(uxy_matrix);

% Get degrees of polynomial d(x,y)
t1 = m1-(m1_k1);
t2 = m2-(m2_k2);

% %
% %
% Build matrices C(u) and C(v)
% Get number of zeros in d(x,y)
nNonZeros_dxy = GetNumNonZeros(t1,t2,k);
nZeros_dxy = (t1+1) * (t2+1) - nNonZeros_dxy;

nNonZeros_ud = GetNumNonZeros(m1,m2,m);
nNonZeros_vd = GetNumNonZeros(n1,n2,n);

C1 = BuildT1_Both(uxy_matrix,m-k,k,t1,t2);
C2 = BuildT1_Both(vxy_matrix,n-k,k,t1,t2);

C = [C1 ; C2];

% % 
% % 
% Preprocess f(x,y) and g(x,y) and get in vector form

% Get fww_matrix as a vector
fxy_vec = GetAsVector(fxy_matrix);

% Remove the zeros from f(x,y) 
fxy_vec = fxy_vec(1:nNonZeros_ud,:);

% get gww_matrix as a vector
gxy_vec = GetAsVector(gxy_matrix);

% Remove the zeros from g(x,y)
gxy_vec = gxy_vec(1:nNonZeros_vd,:);

% Build the RHS vector
rhs_vec = [fxy_vec;
           gxy_vec];

%% Calculate the solution vector

% Calculate the x vector by pinv       
x = SolveAx_b(C,rhs_vec);


% Append the removed zeros
dxy_vec = ...
    [
        x;
        zeros(nZeros_dxy,1);
        ];

% Arrange d(x,y) into a matrix form based on its dimensions.
dxy_matrix = GetAsMatrix(dxy_vec,t1,t2);



end