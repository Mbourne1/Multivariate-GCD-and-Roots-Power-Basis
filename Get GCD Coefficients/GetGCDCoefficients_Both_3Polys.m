function dxy_matrix = GetGCDCoefficients_Both_3Polys(fxy, gxy, hxy, uxy, vxy, wxy, m, n, o, k)
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
[m1, m2] = GetDegree(fxy);

% Get degree of polynomial g(x,y)
[n1, n2] = GetDegree(gxy);

% Get degree of polynomial h(x,y)
[o1, o2] = GetDegree(hxy);

% Get degrees of polynomial u(x,y)
[m1_k1,m2_k2] = GetDegree(uxy);

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
nNonZeros_wd = GetNumNonZeros(o1,o2,o);

C1 = BuildT1_Both(uxy,m-k,k,t1,t2);
C2 = BuildT1_Both(vxy,n-k,k,t1,t2);
C3 = BuildT1_Both(wxy,o-k,k,t1,t2);

C = [C1; C2; C3];

% % 
% % 
% Preprocess f(x,y) and g(x,y) and get in vector form

% Get fww_matrix as a vector
fxy_vec = GetAsVector(fxy);

% Remove the zeros from f(x,y) 
fxy_vec = fxy_vec(1:nNonZeros_ud,:);

% Get gww_matrix as a vector
gxy_vec = GetAsVector(gxy);

% Remove the zeros from g(x,y)
gxy_vec = gxy_vec(1:nNonZeros_vd,:);

% Get hww_matrix as a vector
hxy_vec = GetAsVector(hxy);
hxy_vec = hxy_vec(1:nNonZeros_wd,:);

% Build the RHS vector
rhs_vec = [fxy_vec;
           gxy_vec;
           hxy_vec];

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