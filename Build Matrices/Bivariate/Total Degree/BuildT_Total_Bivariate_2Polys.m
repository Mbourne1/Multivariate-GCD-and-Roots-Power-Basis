function Sk = BuildT_Total_Bivariate_2Polys(fxy,gxy,m,n,k)
% Build the kth Sylvester matrix where polynomials f(x,y) and g(x,y) are
% given in terms of their total degree.
%
%
% % Inputs
%
% [fxy, gxy] : Coefficients of polynomial f(x,y) and g(x,y)
%
% [m, n] : Total degree of f(x,y) and g(x,y)
%
% k : Total degree of d(x,y) and index of the kth Sylvester subresultant
% matrix to be built.
%
% % Outputs.
%
% Sk : The kth Sylvester matrix S_{k}

% Build the partition T_{n-k}(f)
Tf = BuildT1_Total_Bivariate(fxy, m, n-k);

% Build the partiton T_{m-k}(g)
Tg = BuildT1_Total_Bivariate(gxy, n, m-k);

% Build the kth Sylvester subresultant matrix.
Sk = [Tf Tg]; 

end