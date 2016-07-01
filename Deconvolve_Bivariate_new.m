function [hxy_matrix] = Deconvolve_Bivariate_new(fxy_matrix,gxy_matrix,m,n)
% Return the matrix of coefficients of the polynomial h, where h = f/g

%

% Get degrees of polynomial f
[m1,m2] = GetDegree(fxy_matrix);

% Get the degrees of polynomial g
[n1,n2] = GetDegree(gxy_matrix);

% h has total degree m-n
degt_hxy = m - n;
deg1_hxy = m1 - n1;
deg2_hxy = m2 - n2;

nNonZeros_hxy = GetNumNonZeros(deg1_hxy,deg2_hxy,degt_hxy);
nZeros_hxy = (m1-n1 + 1) * (m2-n2 + 1) - nNonZeros_hxy;

% Build the matrix C(g)
C1 = BuildT1(gxy_matrix,m1-n1,m2-n2);

% Remove columns of C(g) corresponding to zeros in the vector of h(x)
C1 = C1(:,1:nNonZeros_hxy);

% Get the polynomial f in vector form
f = GetAsVector(fxy_matrix);

% Solve the Ax=b problem
x_ls = SolveAx_b(C1,f);

hxy_vec = ...
    [
    x_ls;
    zeros(nZeros_hxy,1);
    ];

% Get h(x,y) as a vector
hxy_matrix = GetAsMatrix(hxy_vec,m1-n1,m2-n2);


end

