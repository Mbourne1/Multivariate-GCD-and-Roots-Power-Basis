function Sk = BuildSylvesterMatrix_Total_3Polys(fxy, gxy, hxy, m, n, o, k)
% Build the kth Sylvester matrix where polynomials f(x,y) and g(x,y) are
% given in terms of their total degree.
%
% [fxy, gxy, hxy] : Coefficients of polynomial f(x,y) g(x,y) and h(x,y)
%
% [m, n, o] : Total degree of f(x,y) g(x,y) and h(x,y)
%
% k : Total degree of d(x,y) and index of the kth Sylvester subresultant
% matrix to be built.
%
% % Outputs.
%
% Sk : The kth Sylvester matrix S_{k}

% Build the partition T_{n-k}(f)
T1 = BuildT1_Total(fxy,m,n-k);

T2 = BuildT1_Total(fxy,m,o-k);

% Build the partiton T_{m-k}(g)
T3 = BuildT1_Total(gxy,n,m-k);

T4 = BuildT1_Total(hxy,o,m-k);

diagonal = blkdiag(T1,T2);
column = [T3 ; T4];

% Build the kth Sylvester subresultant matrix.
Sk = [diagonal column]; 

end