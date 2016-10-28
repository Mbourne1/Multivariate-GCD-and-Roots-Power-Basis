function Y = BuildY(m1,m2,n1,n2,t1,t2,opt_col,x_ls,alpha,theta1,theta2)
% This function builds the matrix Y_{t_{1},t_{2}}, 
% Where Y(x)*z = E_{t_{1},t_{2}}(z)*x

% Insert a zero into the position of the optimal_column
first_part = x_ls(1:(opt_col-1));
second_part = x_ls(opt_col:end);
x = [first_part ; 0 ; second_part];

% Separate the x into coefficients of u and coefficients of v

% The first (n1-t1+1) x (n2-t2+1) coefficients are of v
num_coeff_v = (n1-t1+1) * (n2-t2+1);
num_coeff_u = (m1-t1+1) * (m2-t2+1);

% Get the vector of coefficients of v(w,w)
vww_vec = x(1:num_coeff_v);

% Get the vector of coefficients of u(w,w)
uww_vec = - x(num_coeff_v+1:end);


vww_mat = GetAsMatrix(vww_vec,n1-t1,n2-t2);
uww_mat = GetAsMatrix(uww_vec,m1-t1,m2-t2);

T1 = BuildT1(vww_mat,m1,m2);
T2 = BuildT1(-uww_mat,n1,n2);

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