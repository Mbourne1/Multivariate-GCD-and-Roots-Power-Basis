function Yt = BuildY_BothDegree_STLN(x,m,m1,m2,n,n1,n2,k,k1,k2)
% Build the matrix Y_{t}
% Where Y(x) * z = E(z) * x
% The vector x only contains the non-zero entries of v(x,y) and u(x,y)
%
% Inputs
%
% x
%
% m : Total degree of polynomial f(x,y)
%
% m1 : Degree of polynomial f(x,y) with respect to x
%
% m2 : Degree of polynomial f(x,y) with respect to y
%
% n : Total degree of polynomial g(x,y)
%
% n1 : Degree of polynomial g(x,y) with respect to x
%
% n2 : Degree of polynomial g(x,y) with respect to y
%
% k : Total degree of polynomial d(x,y)
%
% k1 : Degree of polynomial d(x,y) with respect to x
%
% k2 : Degree of polynomial d(x,y) with respect to y


% Get number of coefficients in u(x,y)
nCoefficients_uxy = (m1-k1+1) * (m2-k2+1);

% Get number of coefficients in v(x,y)
nCoefficients_vxy = (n1-k1+1) * (n2-k2+1);

% Get number of non-zero Coefficients in u(x,y)
nNonZeros_uxy = GetNumNonZeros(m1-k1,m2-k2,m-k);

% Get number of non-zero coefficients in v(x,y)
nNonZeros_vxy = GetNumNonZeros(n1-k1,n2-k2,n-k);

% Get number of zero coefficients in u(x,y)
nZeros_uxy = nCoefficients_uxy - nNonZeros_uxy;

% Get number of zero coefficients in v(x,y)
nZeros_vxy = nCoefficients_vxy - nNonZeros_vxy;

% Get number of non-zero coefficients in f(x,y)
nNonZeros_fxy = GetNumNonZeros(m1,m2,m);

% Get number of non-zero coefficients in g(x,y)
nNonZeros_gxy = GetNumNonZeros(n1,n2,n);

% % Split the vector x
xv = x(1:nNonZeros_vxy);
xu = 1.* x(nNonZeros_vxy + 1 : end);

% get x_u as a matrix
mat_xu = GetAsMatrix(...
    [
        xu;
        zeros(nZeros_uxy,1)
    ]...
    ,m1-k1,m2-k2);

% Get x_v as a matrix
mat_xv = GetAsMatrix(...
    [
        xv;
        zeros(nZeros_vxy,1);
    ]...
    ,n1-k1,n2-k2);

% Build the matrix C1 and C2
Tv = BuildT1(mat_xv,m1,m2);
Tu = BuildT1(mat_xu,n1,n2);

% Remove the columns of C1(v) and C2(u) corresponding to zeros in f(x,y) and
% g(x,y)
Tv = Tv(:,1:nNonZeros_fxy);
Tu = Tu(:,1:nNonZeros_gxy);

% Remove the rows of C1(v) and C2(u) corresponding to zeros in the product
% f(x,y)*v(x,y) and g(x,y)*u(x,y)

nNonZeros_fv = GetNumNonZeros(m1+n1-k1,m2+n2-k2,m+n-k);
nNonZeros_gu = GetNumNonZeros(n1+m1-k1,n2+m2-k2,n+m-k);

Tv = Tv(1:nNonZeros_fv,:);
Tu = Tu(1:nNonZeros_gu,:);

Yt = [Tv Tu];



end