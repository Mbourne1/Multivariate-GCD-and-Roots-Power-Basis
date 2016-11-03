function Skk1k2 = BuildSylvesterMatrix_Both(fxy,gxy,m,n,k,k1,k2)
% Build the kth Sylvester subresultant matrix S_{k,k1,k2}
%
% % Inputs
%
% fxy : Coefficients of polynomial f(x,y)
%
% gxy : Coefficients of polynomial g(x,y)
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
[m1,m2] = GetDegree(fxy);

% Get the degree of polynomial g(x,y).
[n1,n2] = GetDegree(gxy);

% % Build the partitions of the Sylvester matrix S_{t}

T1 = BuildT1_Both(fxy,m,n-k,n1-k1,n2-k2);

T2 = BuildT1_Both(gxy,n,m-k,m1-k1,m2-k2);

Skk1k2 = [T1 T2];


end
