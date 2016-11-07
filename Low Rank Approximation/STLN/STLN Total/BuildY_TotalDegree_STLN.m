function Y = BuildY_TotalDegree_STLN(x,m,n,k)
% Build the matrix Y_{t}
% Where Y(x) * z = E(z) * x
%
% % Inputs.
%
% x : Vector x least squares solution with zero inserted.
%
% m : Degree of f(x,y)
% 
% n : Degree of g(x,y)
%
% k : Degree of d(x,y)


% Get the  number of coefficients of x_{1}
nCoeff_x1 = nchoosek(n-k+2,2);
nCoeff_x2 = nchoosek(m-k+2,2);

nZeros_x1 = nchoosek(n-k+1,2);
nZeros_x2 = nchoosek(m-k+1,2);

% Get vector of coefficients of x_{1}(x,y)
x1 = x(1:nCoeff_x1);

% Get vector of coefficients of x_{2}(x,y)
x2 = x(nCoeff_x1+1:nCoeff_x1 + nCoeff_x2);

% Get x_{u}(x,y) and x_{v}(x,y) as a matrix
mat_x1 = GetAsMatrix([x1; zeros(nZeros_x1,1)],n-k,n-k);
mat_x2 = GetAsMatrix([x2; zeros(nZeros_x2,1)],m-k,m-k);

% Build the matrices C(v) and C(u)
C1 = BuildT1_Total(mat_x1,n-k,m);
C2 = BuildT1_Total(mat_x2,m-k,n);

% Build the Matrix

Y = [C1 C2];



end