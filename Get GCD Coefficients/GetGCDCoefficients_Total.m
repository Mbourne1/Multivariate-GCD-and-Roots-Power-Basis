function dxy_matrix = GetGCDCoefficients_Total(fxy_matrix,gxy_matrix,uxy_matrix,vxy_matrix,m,n,t)
% Given the matrices of coefficients of f(x,y) and g(x,y), the quotient
% polynomials u(x,y) and v(x,y), and optimal values for alpha, theta_{1}
% and theta_{2}, calculate the coefficients of the GCD d(x,y).
%
% % Inputs
%
% fxy : Coefficients of polynomial f(x,y)
%
% gxy : Coefficients of polynomial g(x,y)
%
% uxy : Coefficients of polynomial u(x,y)
%
% vxy : Coefficients of polynomial v(x,y)
% 
% m : Total degree of polynomial f(x,y)
%
% n : Total degree of polynomial g(x,y)
%
% t : Total degree of polynomial d(x,y)


% % Build Matrix C
% Build the Cauchy matrix of coefficients of u(w,w)
C1 = BuildT1_Total(uxy_matrix,m-t,t);

% Build the Cauchy matrix of coefficients of v(w,w)
C2 = BuildT1_Total(vxy_matrix,n-t,t);

% Build the RHS vector of coefficients of f and g
C = [C1;C2];

% % Build Vector f(w,w)
% Get f with thetas included
% padd f so that it is in terms of its total degree
[m1,m2] = GetDegree(fxy_matrix);
padd = zeros(m+1,m+1);
padd(1:m1+1,1:m2+1) = fxy_matrix;
fxy_matrix = padd;

% % Build Vector g(w,w)
% Get g with thetas included
% padd g so that it is in terms of its total degree
[n1,n2] = GetDegree(gxy_matrix);
padd = zeros(n+1,n+1);
padd(1:n1+1,1:n2+1) = gxy_matrix;
gxy_matrix = padd;

% Get fxy_matrix as a vector
fxy_vec = GetAsVector(fxy_matrix);

% Remove the zeros associated with the polynomial by total degree
num_elements_f = nchoosek(m+2,2);
fxy_vec = fxy_vec(1:num_elements_f);

% get gww_matrix as a vector
gxy_vec = GetAsVector(gxy_matrix);

% Remove the zeros associated with the polynomial by total degree
num_elements_g = nchoosek(n+2,2);
gxy_vec = gxy_vec(1:num_elements_g);

% % Build the RHS vector
rhs_vec = [fxy_vec;
           gxy_vec];

% % Calculate the solution vector

% Calculate the x vector by pinv 
x = SolveAx_b(C,rhs_vec);
    
dxy_vec = x;

% Get the residual associated with the solution x. (Small residual implies good approximation)    
residual = pinv(C)*rhs_vec - x;

% Padd d(w,w) with zeros so that it can be put back into matrix form
dxy_vec = [dxy_vec ; zeros(nchoosek(t+2-1,2),1)];

% Arrange dw into a matrix form based on its dimensions.
dxy_matrix = GetAsMatrix(dxy_vec,t,t);


end