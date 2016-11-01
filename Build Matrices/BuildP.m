function P = BuildP(m1,m2,n1,n2,alpha,theta1,theta2,opt_col,t1,t2,num_cols_T1)
% Calculate the matrix DP where P is the matrix such that c = P[f;g]

% Get the number of coefficients in polynomial f
num_coeff_f = (m1+1).*(m2+1);

% Get the number of coefficients in polynomial g
num_coeff_g = (n1+1).*(n2+1);


    if opt_col <= num_cols_T1
        % Optimal column in first partition
        
        % % Build the matrix P
        
        % Build the matrix P1
        P1 = BuildP_sub(m1,m2,n1,n2,theta1,theta2,opt_col,t1,t2);
        
        % Build the matrix P2
        rows = (m1+n1-t1+1)*(m2+n2-t2+1);
        P2 = zeros(rows,num_coeff_g);
        
        % Build the matrix P
        P = [P1 P2];
        
    else
        % Optimal column in second partition
        
        % Build the matrix P1
        rows = (m1+n1-t1+1)*(m2+n2-t2+1);
        P1 = zeros(rows,num_coeff_f);
        
        % Build the matrix P2
        % Get the position of the optimal column with respect to T(g)
        opt_col_rel = opt_col - num_cols_T1;
        P2 = BuildP_sub(n1,n2,m1,m2,theta1,theta2,opt_col_rel,t1,t2);
        
        % Build the matrix P.
        P = [P1 alpha.*P2];
        
    end
    
end
    
function P = BuildP_sub(m1,m2,n1,n2,theta1,theta2,opt_col,t1,t2)
% Build the matrix P, used in SNTLN function. P is a matrix which is
% obtained from the decomposition of a column vector c_{t} into a matrix 
% vector product P_{t} [f;g]
% c_{t} is the column of the Sylvester matrix S_{t}(f,g).
%
%   %   Inputs
%
%   m1 :    Degree of polynomial f(x,y) with respect to x
%
%   m2 :    Degree of polynomial f(x,y) with respect to y
%
%   n1 :    Degree of polynomial g(x,y) with respect to x
%
%   n2 :    Degree of polynomial g(x,y) with respect to y
%
%   theta1 : Optimal value of theta_{1}
%
%   theta2 : Optimal value of theta_{2}
%
%   opt_col : Optimal column for removal from S(f,g)
%
%   t1 : Degree of GCD d(x,y) with respect to x
%
%   t2 : Degree of GCD d(x,y) with respect to y


% Build the coefficient matrix of thetas 
pre_thetas = diag(theta1.^(0:1:m1));
post_thetas = diag(theta2.^(0:1:m2));
mat = ones(m1+1,m2+1);
mat = pre_thetas * mat * post_thetas;

% pad with zeros
num_mult_wrt_x = n1-t1;
num_mult_wrt_y = n2-t2;

% Produce a zero matrix to fill the space
padd_mat = zeros(m1+n1-t1+1, m2+n2-t2+1);

% % from the index of the optimal column, Get the number of multiplications
% with respec to x and number with respect to y.
[i,j] = GivenCol_getIndex2(n1-t1,n2-t2,opt_col);


ihat = i+1;
jhat = j+1;

num_rows_f = m1+1;
num_cols_f = m2+1;

% inser the theta matrix into the zero matrix
padd_mat(ihat:i+num_rows_f, jhat:j+num_cols_f) = mat;

% Get the padded matrix as a vector
vec_padd_mat = GetAsVector(padd_mat);

% Diagonalise the vector.
diag_mat_vec_padd_mat = diag(vec_padd_mat);


% Remove the zero columns
diag_mat_vec_padd_mat( :, ~any(diag_mat_vec_padd_mat,1) ) = [];  %columns


P = diag_mat_vec_padd_mat;


end