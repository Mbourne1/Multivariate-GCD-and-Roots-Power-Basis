function Sk = BuildT_Relative_Bivariate_3Polys(fxy, gxy, hxy, k1,k2)
% Given two input polynomials f(x,y) and g(x,y), build the (k1,k2)-th
% Sylvester subresultant.
%
% Inputs
%
% [fxy, gxy, hxy] : Coefficients of input polynomial f(x,y), g(x,y) and h(x,y)
%
% k1 : With respect to x.
%
% k2 : With respect to y.
%
% Outputs.
%
% Sk : The Sylvester Subresultant S_{k_{1},k_{2}}


% Get degrees m1 and m2 of polynomial f(x,y)
[m1, m2] = GetDegree_Bivariate(fxy);

% Get degrees n1 and n2 of polynomial g(x,y)
[n1, n2] = GetDegree_Bivariate(gxy);

% Get degrees o1 and o2 of polynomial h(x,y)
[o1, o2] = GetDegree_Bivariate(hxy);

% Build the partitions of the Sylvester subresultant matrix S_{k1,k2}(f,g)
T1 = BuildT1_Relative_Bivariate(fxy, n1-k1, n2-k2);

T2 = BuildT1_Relative_Bivariate(fxy, o1-k1, o2-k2);

T3 = BuildT1_Relative_Bivariate(gxy, m1-k1, m2-k2);
T4 = BuildT1_Relative_Bivariate(hxy, m1-k1, m2-k2);

diagonal = blkdiag(T1,T2);
column = [T3; T4];

% Build the sylvester matrix
Sk = [diagonal column];


end