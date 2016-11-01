
function Yt = BuildY_RelativeDegree_STLN(x,m1,m2,n1,n2,t1,t2)
% Build the matrix Y_{t}
% Where Y(x) * z = E(z) * x

% Get the  number of coefficients of x_{v}, the perturbations of v(x,y)
nCoefficients_xv = (n1-t1+1) * (n2-t2+1);

% Get vector of coefficients of x_{v}(x,y)
xv = x(1:nCoefficients_xv);
% Get vector of coefficients of x_{u}(x,y)
xu = x(nCoefficients_xv+1:end);

% Get x_{u}(x,y) and x_{v}(x,y) as a matrix
mat_xu = GetAsMatrix(xu,m1-t1,m2-t2);
mat_xv = GetAsMatrix(xv,n1-t1,n2-t2);

% Build the matrices C(v) and C(u)
C1 = BuildT1(mat_xv,m1,m2);
C2 = BuildT1(mat_xu,n1,n2);

% Build the Matrix Y = (C(v) C(u)

Yt = [C1 C2];



end