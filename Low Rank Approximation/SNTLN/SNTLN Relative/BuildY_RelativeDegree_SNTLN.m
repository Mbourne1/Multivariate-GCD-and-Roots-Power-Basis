function Y = BuildY_RelativeDegree_SNTLN(x,m1,m2,n1,n2,k1,k2,alpha,theta1,theta2)
% This function builds the matrix Y_{t_{1},t_{2}}, 
% Where Y(x)*z = E_{t_{1},t_{2}}(z)*x


% Separate the x into coefficients of u and coefficients of v

% The first (n1-t1+1) x (n2-t2+1) coefficients are of v
nCoeff_x1 = (n1-k1+1) * (n2-k2+1);

% Get the vector of coefficients of v(w,w)
x1_vec = x(1:nCoeff_x1);

% Get the vector of coefficients of u(w,w)
x2_vec = x(nCoeff_x1+1:end);


x1_mat = GetAsMatrix(x1_vec,n1-k1,n2-k2);
x2_mat = GetAsMatrix(x2_vec,m1-k1,m2-k2);

T1 = BuildT1_Relative(x1_mat,m1,m2);
T2 = BuildT1_Relative(x2_mat,n1,n2);

% The last (m1-t1+1) x (m2-t2+1) coefficients are of u

% Multiply by the thetas of f and g
th1_mat = diag(theta1.^(0:1:m1));
th2_mat = diag(theta2.^(0:1:m2));

fww_thetas_mat = th1_mat * ones(m1+1,m2+1) * th2_mat;
th_mat1 = GetAsVector(fww_thetas_mat);


th1_mat = diag(theta1.^(0:1:n1));
th2_mat = diag(theta2.^(0:1:n2));
gww_thetas_mat = th1_mat * ones(n1+1,n2+1) * th2_mat;
th_mat2 = GetAsVector(gww_thetas_mat);

% multiply

vec = [th_mat1 ; th_mat2];
th_mat = diag(vec);

% multiply by the alpha of g
Y = [T1 alpha*T2] * th_mat;
end