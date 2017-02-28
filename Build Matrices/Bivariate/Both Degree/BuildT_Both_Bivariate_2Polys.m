function Skk1k2 = BuildT_Both_Bivariate_2Polys(fxy, gxy, m, n, k, k1, k2)
% Build the kth Sylvester subresultant matrix S_{k,k1,k2}
%
% % Inputs
%
% [fxy, gxy] : Coefficients of polynomial f(x,y) and g(x,y)
%
% m : Total degree of polynomial f(x,y)
%
% n : Total degree of polynomial g(x,y)
%
% k : Total degree of 
%
% % Outputs
%
% Skk1k2 : Sylvester subresultant matrix S_{k,k1,k2}(f,g)


% Get the degree of polynomial f(x,y).
[m1, m2] = GetDegree_Bivariate(fxy);

% Get the degree of polynomial g(x,y).
[n1, n2] = GetDegree_Bivariate(gxy);

% % Build the partitions of the Sylvester matrix S_{t}% 
T1 = BuildT1_Both_Bivariate(fxy, m, n-k, n1-k1, n2-k2);

T2 = BuildT1_Both_Bivariate(gxy, n, m-k, m1-k1, m2-k2);

Skk1k2 = [T1 T2];


end
