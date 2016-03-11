function [uxy_calc_mtrx,vxy_calc_mtrx] = GetQuotients(fxy_matrix,gxy_matrix,t1,t2,opt_alpha,opt_theta_1,opt_theta_2)
% Given two polynomials and the knowledge of the degree of the GCD. Obtain
% the two quotient polynomials u and v.
%
%       Inputs.
%
%   fxy_matrix :
%
%   gxy_matrix :
%
%   t1 :
%
%   t2 :
%
%   opt_alpha :
%
%   opt_theta_1 :
%   
%   opt_theta_2 :


%%
% Get Data

[r,c] = size(fxy_matrix);
m1 = r-1;
m2 = c-1;

[r,c] = size(gxy_matrix);
n1 = r-1;
n2 = c-1;

%% Preprocess
% Add the preprocessors
% multiply fxy by optimal theta 1
% each row of the matrix is multiplied by theta_{1}^{i} where i is
theta1_mat = diag(opt_theta_1.^(0:1:m1));
theta2_mat = diag(opt_theta_2.^(0:1:m2));
fww_matrix = (theta1_mat * fxy_matrix * theta2_mat);


% Preprocess polynomial gxy
theta1_mat = diag(opt_theta_1.^(0:1:n1));
theta2_mat = diag(opt_theta_2.^(0:1:n2));
gww_matrix = (theta1_mat * gxy_matrix * theta2_mat);


% % Build the partitions of the Sylvester matrix S_{t}

% Build the first partition containing coefficients of fxy
T1 = BuildT1(fww_matrix,n1-t1,n2-t2);

% Build the second partition containing coefficients of gxy
T2 = BuildT1(gww_matrix,m1-t1,m2-t2);

% Concatenate the two partitions
St = [T1 opt_alpha.*T2];

num_zeros_remove_T1 = 0;
num_zeros_remove_T2 = 0;

%%%%%%%%%%%%%%%%%%%%%
%% Get the optimal column for removal
opt_col = GetOptimalColumn(fww_matrix,opt_alpha.*gww_matrix,t1,t2);


%% Having found the optimal column, obtain u and v the quotient polynomials.
Atj = St;
cki = St(:,opt_col);
Atj(:,opt_col) = [];

% Perform QR decomposition of Ak to obtain the solution x
[~,n_col] = size(Atj);
[Q,R] = qr(Atj);
R1 = R(1:n_col,:);
cd = Q'*cki;
c = cd(1:n_col,:);
x_ls = R1\c;



% Obtain the solution vector x = [-v;u]
vecx =[
    x_ls(1:(opt_col)-1);
    -1;
    x_ls(opt_col:end);
    ]  ;

vecx = vecx./vecx(1);

% Get number of coefficients in u(x,y) and v(x,y)
num_coeff_v = (n1-t1+1) * (n2-t2+1);
num_coeff_u = (m1-t1+1) * (m2-t2+1);

% get the vector of coefficients of v
vww_calc = [...
            vecx(1:(num_coeff_v - num_zeros_remove_T1));
            zeros(num_zeros_remove_T1,1)
          ];
      
% get the vector of coefficients of u
uww_calc = [...
            (-1).*vecx((num_coeff_v - num_zeros_remove_T1)+1:end);
            zeros(num_zeros_remove_T2,1);
            ];
        

%% Get u and v in matrix form
% Arrange uw into a matrix form based on their dimensions.

uww_calc_mtrx = GetAsMatrix(uww_calc,m1-t1,m2-t2);
vww_calc_mtrx = GetAsMatrix(vww_calc,n1-t1,n2-t2);

%% Get u(x,y) and v(x,y) from u(w,w) and v(w,w)

th1 = diag(1./opt_theta_1.^(0:1:m1-t1));
th2 = diag(1./opt_theta_2.^(0:1:m2-t2));

uxy_calc_mtrx = th1 * uww_calc_mtrx * th2;

th1 = diag(1./opt_theta_1.^(0:1:n1-t1));
th2 = diag(1./opt_theta_2.^(0:1:n2-t2));

vxy_calc_mtrx = th1 * vww_calc_mtrx * th2;






end
