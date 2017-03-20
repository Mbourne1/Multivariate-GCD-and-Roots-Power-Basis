function [hxy_matrix] = Deconvolve_Bivariate_Single_Respective(fxy, gxy)
% Return the matrix of coefficients of the polynomial h, where h = f/g

%

% Get degrees of polynomial f(x,y)
[m1, m2] = GetDegree_Bivariate(fxy);

% Get the degrees of polynomial g(x,y)
[n1, n2] = GetDegree_Bivariate(gxy);

% Build the matrix C(g)
C1 = BuildT1_Relative_Bivariate(gxy, m1-n1, m2-n2);

% Get the polynomial f(x,y) in vector form
f = GetAsVector(fxy);

% Solve the Ax=b problem
h = SolveAx_b(C1,f);

% Get h(x,y) as a vector
hxy_matrix = GetAsMatrix(h, m1-n1, m2-n2);


end

