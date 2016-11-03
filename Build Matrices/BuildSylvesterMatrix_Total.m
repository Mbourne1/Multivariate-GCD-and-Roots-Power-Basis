function Sk = BuildSylvesterMatrix_Total(fxy,gxy,m,n,k)
% Build the kth Sylvester matrix where polynomials f(x,y) and g(x,y) are
% given in terms of their total degree.
%
% fxy : Coefficients of polynomial f(x,y)
%
% gxy : Coefficients of polynomial g(x,y)
%
% m : Total degree of f(x,y)
%
% n : Total degree of g(x,y)
%
% k : Total degree of d(x,y) and index of the kth Sylvester subresultant
% matrix to be built.
%
% % Outputs.
%
% Sk : The kth Sylvester matrix S_{k}

% Build the partition T_{n-k}(f)
Tf = BuildT1_TotalDegree(fxy,m,n-k);

% Build the partiton T_{m-k}(g)
Tg = BuildT1_TotalDegree(gxy,n,m-k);

% Build the kth Sylvester subresultant matrix.
Sk = [Tf Tg]; 

end