function Skk1k2 = BuildSylvesterMatrix_Both_3Polys(fxy, gxy, hxy, m, n, o, k, k1, k2)
% Build the kth Sylvester subresultant matrix S_{k,k1,k2}
%
% % Inputs
%
% [fxy, gxy, hxy] : Coefficients of polynomial f(x,y), g(x,y) and h(x,y)
%
% [m, n, o] : Total degree of polynomial f(x,y) g(x,y) and h(x,y)
%
% k : Total degree of 
%
% % Outputs
%
% Skk1k2 : Sylvester subresultant matrix S_{k,k1,k2}(f,g)


% Get the degree of polynomial f(x,y).
[m1, m2] = GetDegree(fxy);

% Get the degree of polynomial g(x,y).
[n1, n2] = GetDegree(gxy);

% Get the degree of polynomial h(x,y).
[o1, o2] = GetDegree(hxy);


% % Build the partitions of the Sylvester matrix S_{t}

T1 = BuildT1_Both(fxy, m, n-k, n1-k1, n2-k2);
T2 = BuildT1_Both(fxy, m, o-k, o1-k1, o2-k2);


T3 = BuildT1_Both(gxy, n, m-k, m1-k1, m2-k2);
T4 = BuildT1_Both(hxy, o, m-k, m1-k1, m2-k2);

diagonal = blkdiag(T1,T2);
column = [T3;T4];

Skk1k2 = [diagonal column];


end
