function S = BuildSylvesterMatrix_Respective(fxy_matrix,gxy_matrix,k1,k2)
% Given two input polynomials f(x,y) and g(x,y), build the (k1,k2)-th
% Sylvester subresultant.
%
% Inputs
%
% fxy_matrix : Coefficients of input polynomial f(x,y).
%
% gxy_matrix : Coefficients of input polynomial g(x,y).
%
% k1 : With respect to x.
%
% k2 : With respect to y.
%
% Outputs.
%
% S : The Sylvester Subresultant S_{k_{1},k_{2}}


% Get degrees m1 and m2 of polynomial f(x,y)
[m1,m2] = GetDegree(fxy_matrix);

% Get degrees n1 and n2 of polynomial g(x,y)
[n1,n2] = GetDegree(gxy_matrix);

% Build the partitions of the Sylvester subresultant matrix S_{k1,k2}(f,g)
T1 = BuildT1(fxy_matrix,n1-k1,n2-k2);
T2 = BuildT1(gxy_matrix,m1-k1,m2-k2);

% Build the sylvester matrix
S = [T1 T2]; 
end