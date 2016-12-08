function Sk = BuildSylvesterMatrix_Relative_2Polys(fxy, gxy, k1,k2)
% Given two input polynomials f(x,y) and g(x,y), build the (k1,k2)-th
% Sylvester subresultant.
%
% Inputs
%
% fxy : Coefficients of input polynomial f(x,y).
%
% gxy : Coefficients of input polynomial g(x,y).
%
% k1 : With respect to x.
%
% k2 : With respect to y.
%
% Outputs.
%
% Sk : The Sylvester Subresultant S_{k_{1},k_{2}}


% Get degrees m1 and m2 of polynomial f(x,y)
[m1,m2] = GetDegree(fxy);

% Get degrees n1 and n2 of polynomial g(x,y)
[n1,n2] = GetDegree(gxy);

% Build the partitions of the Sylvester subresultant matrix S_{k1,k2}(f,g)
T1 = BuildT1_Relative(fxy,n1-k1,n2-k2);
T2 = BuildT1_Relative(gxy,m1-k1,m2-k2);

% Build the sylvester matrix
Sk = [T1 T2]; 
end