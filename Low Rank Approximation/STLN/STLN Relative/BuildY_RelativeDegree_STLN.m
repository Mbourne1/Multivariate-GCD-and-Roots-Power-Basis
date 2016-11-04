
function Y = BuildY_RelativeDegree_STLN(x,m1,m2,n1,n2,k1,k2)
% Build the matrix Y_{t}
% Where Y(x) * z = E(z) * x
%
% % Inputs.
%
% x : Vector x least squares solution with zero inserted.
%
% m1 : Degree of f(x,y) with respect to x
% 
% m2 : Degree of f(x,y) with respect to y 
% 
% n1 : Degree of g(x,y) with respect to x
%
% n2 : Degree of g(x,y) with respect to y
%
% k1 : Degree of d(x,y) with respect to x
%
% k2 : Degree of d(x,y) with respect to y


% Get the  number of coefficients of x_{1}, the perturbations of v(x,y)
nCoefficients_x1 = (n1-k1+1) * (n2-k2+1);

% Get vector of coefficients of x_{1}(x,y)
x1 = x(1:nCoefficients_x1);

% Get vector of coefficients of x_{2}(x,y)
x2 = x(nCoefficients_x1+1:end);

% Get x_{u}(x,y) and x_{v}(x,y) as a matrix
mat_x1 = GetAsMatrix(x2,m1-k1,m2-k2);
mat_x2 = GetAsMatrix(x1,n1-k1,n2-k2);

% Build the matrices C(v) and C(u)
C1 = BuildT1_Relative(mat_x2,m1,m2);
C2 = BuildT1_Relative(mat_x1,n1,n2);

% Build the Matrix Y = (C(v) C(u)

Y = [C1 C2];



end