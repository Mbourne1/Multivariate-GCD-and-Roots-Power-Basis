function Sk = BuildSylvesterMatrix_Both(fxy,gxy,m,n,k,k1,k2)
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


% Get the degree of polynomial f(x,y).
[m1,m2] = GetDegree(fxy);

% Get the degree of polynomial g(x,y).
[n1,n2] = GetDegree(gxy);

% % Build the partitions of the Sylvester matrix S_{t}

% Build the first partition containing coefficients of fxy
T1 = BuildT1(fxy,n1-k1,n2-k2);

% Build the second partition containing coefficients of gxy
T2 = BuildT1(gxy,m1-k1,m2-k2);

% Get number of non-zero entries in u(x,y) and v(x,y)
nNoneZeros_uxy = GetNumNonZeros(m1-k1,m2-k2,m-k);
nNoneZeros_vxy = GetNumNonZeros(n1-k1,n2-k2,n-k);

% Get number of zero entries in u(x,y) and v(x,y)
nZeros_uxy = (m1-k1+1) * (m2-k2+1) - nNoneZeros_uxy;
nZeros_vxy = (n1-k1+1) * (n2-k2+1) - nNoneZeros_vxy;

% %
% %
% Remove Columns from S(f,g)

% Remove the zero columns from T_{n1-t1,n2-t2}
T1 = T1(:,1:nNoneZeros_vxy);

% Remove the zero columns from T_{m1-t1,m2-t2}
T2 = T2(:,1:nNoneZeros_uxy);

% %
% % 
% Remove Rows from S(f,g)
nNoneZeros_fv = GetNumNonZeros(m1+n1-k1,m2+n2-k2,m+n-k);
nNoneZeros_gu = GetNumNonZeros(n1+m1-k1,n2+m2-k2,n+m-k);

T1 = T1(1:nNoneZeros_fv,:);
T2 = T2(1:nNoneZeros_gu,:);


% Form the Sylvester matrix.
Sk = [T1 T2];


end
