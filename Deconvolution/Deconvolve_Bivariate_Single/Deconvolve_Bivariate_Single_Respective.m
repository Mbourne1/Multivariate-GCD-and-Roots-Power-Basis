function [hxy_matrix] = Deconvolve_Bivariate_Single_Respective(fxy_matrix,gxy_matrix)
% Return the matrix of coefficients of the polynomial h, where h = f/g

%

% Get degrees of polynomial f(x,y)
[m1,m2] = GetDegree(fxy_matrix);

% Get the degrees of polynomial g(x,y)
[n1,n2] = GetDegree(gxy_matrix);

% Build the matrix C(g)
C1 = BuildT1_Relative(gxy_matrix,m1-n1,m2-n2);

% Get the polynomial f(x,y) in vector form
f = GetAsVector(fxy_matrix);

% Solve the Ax=b problem
h = SolveAx_b(C1,f);

% Get h(x,y) as a vector
hxy_matrix = GetAsMatrix(h,m1-n1,m2-n2);


end

