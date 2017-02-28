function [hxy_matrix] = Deconvolve_Bivariate_Single_Total(fxy, gxy, m, n)
% Return the matrix of coefficients of the polynomial h, where h = f/g

%

% Get degrees of polynomial f(x,y)
[m1, m2] = GetDegree_Bivariate(fxy);

% Make sure that f(x,y) is in a matrix of size (m+1) x (m+1) since
% deconvolution is being computed in terms of total degree.
temp_mat = zeros(m+1, m+1);
temp_mat(1:m1+1, 1:m2+1) = fxy;
fxy = temp_mat;


% Get the degrees of polynomial g(x,y)
[n1, n2] = GetDegree_Bivariate(gxy);

% Make sure that g(x,y) is in a matrix of size (n+1) x (n+1) since
% deconvoltuion is being computed in terms of total degree
temp_mat = zeros(n+1, n+1);
temp_mat(1:n1+1, 1:n2+1) = gxy;
gxy = temp_mat;


% Build the matrix C(g)
C1 = BuildT1_Total(gxy, n, m-n);

% Get the polynomial f(x,y) in vector form
v_fxy = GetAsVector(fxy);

% Get number of coefficients in f(x,y)
nCoefficients_fxy = nchoosek(m+2, 2);

% Remove zeros from f(x,y) vector
v_fxy = v_fxy(1:nCoefficients_fxy);

% Solve the Ax=b problem
v_hxy = SolveAx_b(C1,v_fxy);

nCoefficients_hxy = nchoosek(m-n+2,2);

% Append zeros to the vector v_hxy 
temp_vec = zeros(m-n+1,m-n+1);
temp_vec(1:nCoefficients_hxy) = v_hxy;
v_hxy = temp_vec;


% Get h(x,y) as a vector
hxy_matrix = GetAsMatrix(v_hxy, m-n, m-n);


end

