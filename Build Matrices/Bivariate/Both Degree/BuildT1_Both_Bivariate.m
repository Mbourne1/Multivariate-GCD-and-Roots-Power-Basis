function T1 = BuildT1_Both_Bivariate(fxy, m, n_k, n1_k1, n2_k2)
% Get the matrix T_{}(f) where f(x,y) is defined in terms of both total and
% relative degree.
%
% % Inputs.
%
% fxy : Coefficients of polynomial f(x,y)
%
% m : Total degree of polynomial f(x,y)
%
% n_k : Total degree of v(x,y)
%
% n1_k1 : Degree of v(x,y) with respect to x
%
% n2_k2 : Degree of v(x,y) with respect to y

% Get the degree of polynomial f(x,y).
[m1,m2] = GetDegree_Bivariate(fxy);

% Build the first partition containing coefficients of fxy
T1 = BuildT1_Relative_Bivariate(fxy,n1_k1,n2_k2);

% Get number of nonzeros in v(x,y)
nNoneZeros_vxy = GetNumNonZeros(n1_k1,n2_k2,n_k);

% Remove the zero columns from T_{n1-t1,n2-t2}
T1 = T1(:,1:nNoneZeros_vxy);

% Get number of nonzeros in fv(x,y)
nNoneZeros_fv = GetNumNonZeros(m1+n1_k1,m2+n2_k2,m+n_k);

% Remove rows from T1
T1 = T1(1:nNoneZeros_fv,:);


end