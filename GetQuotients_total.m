function [uxy_calc_mtrx,vxy_calc_mtrx] = ...
    GetQuotients_total(fxy_matrix,gxy_matrix,m,n,t,opt_alpha,opt_theta_1,opt_theta_2)
% Given two polynomials and the knowledge of the degree of the GCD. Obtain
% the two quotient polynomials u and v.
%
%       Inputs.
%
%   fxy_matrix :
%
%   gxy_matrix :
%
%   t :
%
%   opt_alpha :
%
%   opt_theta_1 :
%   
%   opt_theta_2 :

% Replace fxy_matrix with the padded version
[r,c] = size(fxy_matrix);
m1 = r - 1;
m2 = c - 1;
padd_matrix = zeros(m+1,m+1);
padd_matrix(1:m1+1,1:m2+1) = fxy_matrix;
fxy_matrix = padd_matrix;

% Replace gxy_matrix with the padded version
[r,c] = size(gxy_matrix);
n1 = r - 1;
n2 = c - 1;
padd_matrix = zeros(n+1,n+1);
padd_matrix(1:n1+1,1:n2+1) = gxy_matrix;
gxy_matrix = padd_matrix;


%% Preprocess
% Add the preprocessors
% multiply fxy by optimal theta 1
% each row of the matrix is multiplied by theta_{1}^{i} where i is
theta1_mat = diag(opt_theta_1.^(0:1:m));
theta2_mat = diag(opt_theta_2.^(0:1:m));
fww_matrix = (theta1_mat * fxy_matrix * theta2_mat);


% Preprocess polynomial gxy
theta1_mat = diag(opt_theta_1.^(0:1:n));
theta2_mat = diag(opt_theta_2.^(0:1:n));
gww_matrix = (theta1_mat * gxy_matrix * theta2_mat);


% % Build the partitions of the Sylvester matrix S_{t}

% Build the first partition containing coefficients of fxy
T1 = BuildT1_totaldegree(fww_matrix,m,n,t);

% Build the second partition containing coefficients of gxy
T2 = BuildT1_totaldegree(gww_matrix,n,m,t);

% Concatenate the two partitions
St = [T1 opt_alpha.*T2];

num_zeros_remove_T1 = 0;
num_zeros_remove_T2 = 0;

%%%%%%%%%%%%%%%%%%%%%
%% Get the optimal column for removal
opt_col = getOptimalColumn_total(fww_matrix,opt_alpha.*gww_matrix,m,n,t);


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
num_coeff_v = nchoosek(n-t+2,2);
num_coeff_u = nchoosek(m-t+2,2);

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


zeros_vww = zeros(nchoosek(n-t-1+2,2),1);
zeros_uww = zeros(nchoosek(m-t-1+2,2),1);

uww_calc_mtrx = getAsMatrix([uww_calc;zeros_uww],m-t,m-t);
vww_calc_mtrx = getAsMatrix([vww_calc;zeros_vww],n-t,n-t);

%% Get u(x,y) and v(x,y) from u(w,w) and v(w,w)

th1 = diag(1./opt_theta_1.^(0:1:m-t));
th2 = diag(1./opt_theta_2.^(0:1:m-t));

uxy_calc_mtrx = th1 * uww_calc_mtrx * th2;

th1 = diag(1./opt_theta_1.^(0:1:n-t));
th2 = diag(1./opt_theta_2.^(0:1:n-t));

vxy_calc_mtrx = th1 * vww_calc_mtrx * th2;






end
