function [hxy_matrix] = Deconvolve_Bivariate_both(fxy_matrix,gxy_matrix,m,n)
% Perform polynomial deconvolution of the bivariate polynomials f(x,y) and 
% g(x,y) where the total degrees m and n of f and g respectively are known.
% This method is referred to as 'both' since BOTH the total degree and
% relative degrees of the polynomials are utilised.
%
% Inputs.
%
% fxy_matrix : Coefficients of polynomial f(x,y)
%
% gxy_matrix : Coefficients of polynomial g(x,y)
%
% m : Total degree of f(x,y)
%
% n : Total degree of g(x,y)
%
%
% Outputs.
%
% hxy_matrix : Coefficients of polynomial h(x,y)

% Return the matrix of coefficients of the polynomial h, where h = f/g


% Get the degrees of polynomial f(x,y)
[m1,m2] = GetDegree(fxy_matrix);

% Get the degrees of polynomial g(x,y)
[n1,n2] = GetDegree(gxy_matrix);

% Get the degrees of polynomial h(x,y)
degt_hxy = m - n;
deg1_hxy = m1 - n1;
deg2_hxy = m2 - n2;

% Get number of zeros in h(x,y)
nNonZeros_hxy = GetNumNonZeros(deg1_hxy,deg2_hxy,degt_hxy);
nZeros_hxy = (m1-n1 + 1) * (m2-n2 + 1) - nNonZeros_hxy;

% Get number of nonzeros in g*h = f(x,y)
nNonZeros_gh = GetNumNonZeros(m1,m2,m);

% % 
% %
% Build the matrix C(g)
C1 = BuildT1(gxy_matrix,m1-n1,m2-n2);

% Remove columns of C(g) corresponding to zeros in the vector of h(x)
C1 = C1(:,1:nNonZeros_hxy);

% Remove the rows of C(g) corresponding to zeros in the vector f(x)
C1 = C1(1:nNonZeros_gh,:);

% %
% %
% Create Right hand side vector f(x,y)

% Get the polynomial f(x,y) in vector form
f = GetAsVector(fxy_matrix);

% Remove zeros from f(x,y).
f = f(1:nNonZeros_gh,1);

% Solve the Ax=b problem.
x_ls = SolveAx_b(C1,f);

% Get set of coefficients of h(x,y), including zeros to form matrix of
% coefficients
hxy_vec = ...
    [
    x_ls;
    zeros(nZeros_hxy,1);
    ];

% Get h(x,y) as a vector
hxy_matrix = GetAsMatrix(hxy_vec,m1-n1,m2-n2);


end

