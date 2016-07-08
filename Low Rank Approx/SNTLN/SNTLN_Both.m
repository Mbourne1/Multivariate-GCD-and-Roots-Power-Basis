function [ o_fxy, o_gxy, o_alpha,o_th1,o_th2,o_X] = ...
    SNTLN_Respective( fxy_matrix,gxy_matrix, i_alpha, i_th1, i_th2,m,n,t,t1,t2,opt_col)
% Obtain the low rank approximation of the Sylvester matrix D*T_{t}(f,g)*Q =
% S_{t}(f,g)
%
% SNTLN( fxy_matrix,gxy_matrix, i_alpha, i_th1, i_th2,t1,t2,opt_col)
%
% Inputs:
%
%
% fxy_matrix : Coefficients of polynomial f, in standard bernstein basis.
%
% gxy_matrix : Coefficients of polynomial g, in standard bernstein basis.
%
% i_alpha : Initial value of alpha
%
% i_th1 : Initial value of theta1
%
% i_th2 : Initial value of theta2
%
% t1 : Degree of AGCD d(x,y) with respect to x
%
% t2 : Degree of AGCD d(x,y) with repsect to y
%
% opt_col : Optimal column for removal from the sylvester matrix, such that col
%           is the column which is most likely a linear combination of the others.
%
%
% Outputs:
%
%
% o_fxy : Coefficients of f(x,y) on output, in standard bernstein basis,
% including added structured perturbations.
%
% o_gxy : Coefficients of fx on output, in standard bernstein basis,
% including added structured perturbations.
%
% o_alpha :
%
% o_th1 :
%
% o_th2 :
%
% o_X :
%

% Global Inputs

global MAX_ERROR_SNTLN
global MAX_ITERATIONS_SNTLN
global PLOT_GRAPHS

% Set the initial iterations number
ite = 1;

% Set initial values of alpha and theta
th1(ite) = i_th1;
th2(ite) = i_th2;
alpha(ite) = i_alpha;

% Get degree of polynomials f.
[m1,m2] = GetDegree(fxy_matrix);

% Get degree of polynomial g.
[n1,n2] = GetDegree(gxy_matrix);

% Get the number of coefficients in the polynomial f(x,y)
num_coeff_f = (m1+1) * (m2+1);

% Get the number of coefficients in the polynomial g(x,y)
num_coeff_g = (n1+1) * (n2+1);

% Get the number of coefficients in both f(x,y) and g(x,y)
num_coeff = num_coeff_f + num_coeff_g;

% Get the number of coefficients in v(x,y)
num_coeff_v = (n1-t1+1) * (n2-t2+1);

% Get the number of coefficients in u(x,y)
num_coeff_u = (m1-t1+1) * (m2-t2+1);

% Get the number of coefficients in the unknown vector x, where A_{t}x =
% c_{t}.
num_coeff_x = num_coeff_u + num_coeff_v - 1;

% Create the identity matrix I, the matrix M formed from I by removing the
% column equivalent to the optimal column for removal from the Sylvester
% subresultant matrix, so that S_{t}(f,g)*M = A_{t}, where A_{t} is the
% Sylvester subresultant matrix with the column removed.

% Get the number of columns in C_{t}(f), the first partition of the Sylvester
% Matrix S_{t}(f,g)
num_cols_in_f_partition = (n1-t1+1) * (n2-t2+1);

% Get the number of columns in C_{t}(g), the second partition of the
% Sylvester matrix S_{t}(f,g)
num_cols_in_g_partition = (m1-t1+1) * (m2-t2+1);

% Get the total number of columns in the Sylvester matrix S_{t}(f,g)
num_cols = num_cols_in_f_partition + num_cols_in_g_partition;

% Create the identity matrix
I = eye(num_cols, num_cols);

% Create the matrix M, such that S(f,g)*M gives A_{t}, the Sylvester Matrix
% with the optimal column removed.
M = I;
M(:,opt_col) = [];

% Let e be the column removed from the identity matrix, such that
% S_{t}(f,g) * e gives the column c_{t}, where c_{t} is the optimal column
% removed from the Sylvester subresultant.
e = I(:,opt_col);

%% Preprocessing
% Obtain polynomials in Modified Bernstein Basis, using initial values of
% alpha and theta.

% Multiply the rows of fxy_matrix by theta1, and multiply the cols of
% fxy_matrix by theta2
fww_matrix = GetWithThetas(fxy_matrix,th1(ite),th2(ite));

% Multiply the rows of gxy_matrix by theta1, and multiply the cols of
% gxy_matrix by theta2
gww_matrix = GetWithThetas(gxy_matrix,th1(ite),th2(ite));

% Form the Coefficient Matrix T = [C(f(w1,w2))|alpha * C(g(w1,w2))] such that T*x = [col]
T1 = BuildT1(fww_matrix,n1-t1,n2-t2);
T2 = BuildT1(gww_matrix,m1-t1,m2-t2);
T = [T1 alpha(ite).*T2];

%
% Calculate the partial derivatives of f(w,w) and g(w,w) with respect to \alpha
Partial_fw_wrt_alpha            = zeros(m1+1,m2+1);
Partial_alpha_gw_wrt_alpha      = gxy_matrix;

%
% Calculate the partial derivatives of f(w,w) with respect to \theta_1
theta_mat = diag((0:1:m1) ./ th1(ite));
Partial_fw_wrt_theta1    = theta_mat * fww_matrix;

%
% Calculate the partial derivative of g(w,w) with respect to theta_1
theta_mat = diag((0:1:n1) ./ th1(ite));
Partial_gw_wrt_theta1    = theta_mat * gww_matrix;

%
% Calculate the partial derivative of f(w,w) with respect to theta_2
theta_mat = diag((0:1:m2) ./ th2(ite));
Partial_fw_wrt_theta2 = fww_matrix * theta_mat;

%
% Calculate the partial deriviates of g(w,w) with respect to theta_2
theta_mat = diag((0:1:n2)./ th2(ite));
Partial_gw_wrt_theta2 = gww_matrix * theta_mat;

%%
% Build the derivative of T(f,g) with respect to alpha
T1_wrt_alpha = BuildT1(Partial_fw_wrt_alpha,n1-t1,n2-t2);
T2_wrt_alpha = BuildT1(Partial_alpha_gw_wrt_alpha,m1-t1,m2-t2);
T_alpha = [T1_wrt_alpha T2_wrt_alpha];

%%
% Calculate the derivative of T(f,g) with respect to theta_{1}
T1_wrt_theta1 = BuildT1(Partial_fw_wrt_theta1,n1-t1,n2-t2);
T2_wrt_theta1 = BuildT1(Partial_gw_wrt_theta1,m1-t1,m2-t2);
Partial_T_wrt_theta1 = [T1_wrt_theta1 alpha(ite)* T2_wrt_theta1];

%%
% Calcualte the derivative of T(f,g) with respect to theta_2
T1_wrt_theta2 = BuildT1(Partial_fw_wrt_theta2,n1-t1,n2-t2);
T2_wrt_theta2 = BuildT1(Partial_gw_wrt_theta2,m1-t1,m2-t2);
Partial_T_wrt_theta2 = [T1_wrt_theta2 alpha(ite)*T2_wrt_theta2];

%%
% Initialise the vector z of structured perturbations
% if we are working with strictly the roots problem, the number of entries
% in z can be reduced.

num_rows_Sylv_mat = (m1 + n1 - t1 + 1) * (m2 + n2 - t2 + 1);
num_cols_T1 = (n1 - t1 + 1) * (n2 - t2 + 1);
num_cols_T2 = (m1 - t1 + 1) * (m2 - t2 + 1);

zk = zeros(num_coeff , 1);

%%
% Initilaise the derivative of N wrt alpha.
Partial_N_wrt_alpha   = zeros(num_rows_Sylv_mat,num_cols_T1 + num_cols_T2);

% Initilaise the derivative of N wrt theta_1.
Partial_N_wrt_theta1   = zeros(num_rows_Sylv_mat,num_cols_T1 + num_cols_T2);

% Initialise the derivative of N wrt theta 2
Partial_N_wrt_theta2   = zeros(num_rows_Sylv_mat,num_cols_T1 + num_cols_T2);

%%
% Initialise the derivative of h
% Calculate the derivatives wrt alpha and theta of the column of DNQ
% that is moved to the right hand side.
Partial_h_wrt_alpha     = Partial_N_wrt_alpha*e;
Partial_h_wrt_theta1    = Partial_N_wrt_theta1*e;
Partial_h_wrt_theta2    = Partial_N_wrt_theta2*e;

% Get the coefficients of f(x,y) in matrix form
fxy_vec = GetAsVector(fxy_matrix);
gxy_vec = GetAsVector(gxy_matrix);

% Get the matrix A_{k}(f,g), which is the subresultant matrix S(f,g) with
% an opitmal column removed
Ak = T;
ck = T(:,opt_col);
Ak(:,opt_col) = [];

%% Build the matrix P
P = BuildP(m1,m2,n1,n2,alpha(ite),th1(ite),th2(ite),opt_col,t1,t2,num_cols_T1);
% Test P
%ck;
%ck2 = P*[fxy_vec;gxy_vec];
%ck - ck2;

%%
% Calculate the derivatives of ck wrt alpha and theta.
Partial_ck_wrt_alpha        = T_alpha*e;
Partial_ck_wrt_theta1       = Partial_T_wrt_theta1*e;
Partial_ck_wrt_theta2       = Partial_T_wrt_theta2*e;


%%
% Perform QR decomposition of Ak to obtain the solution x
x_ls = SolveAx_b(Ak,ck);

% store this version of x_ls
initial_xls = x_ls;
first_part = x_ls(1:(opt_col-1));
second_part = x_ls(opt_col:end);
x = [first_part ; 0 ; second_part];

% Build Matrix Y, where Y(v,u)*[f;g] = S(f,g)*[u;v]
Y = BuildY(m1,m2,n1,n2,t1,t2,opt_col,x_ls,alpha(ite),th1(ite),th2(ite));
% Test Y
%test1 = Y * [fxy_vec;gxy_vec];
%test2 = T * x;
%test1-test2;

% Calculate the initial residual r = ck - (Ak*x)
res_vec = ck - (T*M*x_ls);

% Get the matrix p, which will store all the perturbations returned from LSE file
num_entries = num_coeff_f...
    + num_coeff_g ...
    + num_coeff_x ...
    + 3;


% Set the initial value of vector p to be zero
f = zeros(num_coeff_f...
    + num_coeff_g ...
    + num_coeff_x ...
    + 3,1);

%
% Set the intial value of E to the identity matrix
E = eye(num_entries);

%
% Create the matrix (T+N), initially N is empty so this is the same as T.
T1 = BuildT1(fww_matrix,n1-t1,n2-t2);
T2 = BuildT1(gww_matrix,m1-t1,m2-t2);
TN = [T1 alpha(ite).*T2];

%
% Create The matrix (T+N) with respect to alpha
T1_wrt_alpha = BuildT1(Partial_fw_wrt_alpha,n1-t1,n2-t2);
T2_wrt_alpha = BuildT1(Partial_alpha_gw_wrt_alpha,m1-t1,m2-t2);
TN_wrt_alpha = [T1_wrt_alpha T2_wrt_alpha];

%
% Create The matrix (T+N) with respect to theta1
T1_wrt_theta1 = BuildT1(Partial_fw_wrt_theta1,n1-t1,n2-t2);
T2_wrt_theta1 = BuildT1(Partial_gw_wrt_theta1,m1-t1,m2-t2);
TN_wrt_theta1 = [T1_wrt_theta1 alpha(ite) * T2_wrt_theta1];

%
% Create The matrix (T+N) with respect to theta2
T1_wrt_theta2 = BuildT1(Partial_fw_wrt_theta2,n1-t1,n2-t2);
T2_wrt_theta2 = BuildT1(Partial_gw_wrt_theta2,m1-t1,m2-t2);
TN_wrt_theta2 = [T1_wrt_theta2 alpha(ite) * T2_wrt_theta2];

%%
% Create the matrix C for input into iteration

H_z     = Y-P;

H_x     = TN*M;

H_alpha  = TN_wrt_alpha*M*x_ls - ...
    (Partial_ck_wrt_alpha + Partial_h_wrt_alpha);

H_theta1 = TN_wrt_theta1*M*x_ls - ...
    (Partial_ck_wrt_theta1 + Partial_h_wrt_theta1);

H_theta2 = TN_wrt_theta2*M*x_ls - ...
    (Partial_ck_wrt_theta2 + Partial_h_wrt_theta2);

C       = [H_z H_x H_alpha H_theta1 H_theta2];

% Define the starting vector for the iterations for the LSE problem.
start_point     =   ...
    [...
    zk;...
    x_ls;...
    alpha(ite);...
    th1(ite);...
    th2(ite)
    ];

yy              =   start_point;

% Set the termination criterion to a large value. It will be
% over written later.
condition(ite) = norm(res_vec) ./ norm(ck);

xk = x_ls;


while condition(ite) >(MAX_ERROR_SNTLN) &&  ite < MAX_ITERATIONS_SNTLN
    %while   ite < max_iterations
    
    
    % Use the QR decomposition to solve the LSE problem
    % min |y-p| subject to Cy=q
    y = LSE(E,f,C,res_vec);
    
    % Increment the iteration number
    ite = ite + 1;
    
    % Add the small changes found in LSE problem to existing values
    yy = yy + y;
    
    % obtain the small changes
    
    %% Break down y into its sections
    
    % get the coefficients corresponding to f and g
    delta_zk        = y(1:num_coeff_f + num_coeff_g ,1);
    
    % Remove the zk coefficients from the list of coefficients
    y(1:num_coeff_f + num_coeff_g) = [];
    
    % Get the coefficients corresponding to x
    delta_xk        = y(1:num_coeff_x,1);
    
    % Remove them from the list of coefficients
    y(1:num_coeff_x) = [];
    
    % Get the coefficient corresponding to alpha
    delta_alpha     = y(1:1);
    y(1) = [];
    
    % Get the coefficient corresponding to theta1
    delta_theta1    = y(1:1);
    y(1) = [];
    
    % Get the coefficient corresponding to theta2
    delta_theta2    = y(1:1);
    y(1) = [];
    
    %% Update the variables
    
    % Update variables z_{k}, where z_{k} are perturbations in the
    % coefficients of f and g.
    zk = zk + delta_zk;
    
    % Update x_{k}, where x_{k} is the solution vector, containing
    % coefficients u and v.
    xk = xk + delta_xk;
    
    % Update alpha
    alpha(ite) = alpha(ite-1) + delta_alpha;
    
    % Update theta_{1}
    th1(ite) = th1(ite-1) + delta_theta1;
    
    % Update theta_{2}
    th2(ite) = th2(ite-1) + delta_theta2;
    
    % Obtain polynomials in modified bersntein basis a_{i}\theta^{i}
    
    % Obtain new f(w,w) with improved theta1, and theta2
    fww_matrix = GetWithThetas(fxy_matrix,th1(ite),th2(ite));
    
    % Obtain new g(w,w) with improved theta1 and theta2
    gww_matrix = GetWithThetas(gxy_matrix,th1(ite),th2(ite));
    
    % Construct the Sylvester subresultant matrix S.
    T1 = BuildT1(fww_matrix,n1-t1,n2-t2);
    T2 = BuildT1(gww_matrix,m1-t1,m2-t2);
    T = [T1 alpha(ite).*T2];
    
    % Calculate the partial derivatives of fw and gw with respect to alpha
    Partial_fw_wrt_alpha            = zeros(m1+1,m2+1);
    Partial_alpha_gw_wrt_alpha      = gww_matrix;
    
    %%
    % Calculate the partial derivatives of fw and gw with respect to theta1
    % divide the rows by theta1 and multiply by the old power
    
    % Get the partial derivative of f with respect to theta 1
    temp_mat = diag((0:1:m1)./th1(ite));
    Partial_fw_wrt_theta1 = temp_mat * fww_matrix;
    
    % Get the partial derivative of g with respect to theta1
    temp_mat = diag((0:1:n1)./th1(ite));
    Partial_gw_wrt_theta1 = temp_mat * gww_matrix;
    
    % Get the partial derivative of f with respect to theta2
    temp_mat = diag((0:1:m2)./th2(ite));
    Partial_fw_wrt_theta2 =  fww_matrix * temp_mat;
    
    % Get the partial derivative of g with respect too theta2
    temp_mat = diag((0:1:n2)./th2(ite));
    Partial_gw_wrt_theta2 =  gww_matrix * temp_mat;
    
    
    % Calculate the Partial derivative of T with respect to alpha.
    T1_wrt_alpha = BuildT1(Partial_fw_wrt_alpha,n1-t1,n2-t2);
    T2_wrt_alpha = BuildT1(Partial_alpha_gw_wrt_alpha,m1-t1,m2-t2);
    Partial_T_wrt_alpha = [T1_wrt_alpha T2_wrt_alpha];
    
    
    % Calculate the partial derivative of T with respect to theta1
    T1_wrt_theta1 = BuildT1(Partial_fw_wrt_theta1,n1-t1,n2-t2);
    T2_wrt_theta1 = BuildT1(Partial_gw_wrt_theta1,m1-t1,m2-t2);
    Partial_T_wrt_theta1 = [T1_wrt_theta1 alpha(ite)*T2_wrt_theta1];
    
    % Calculate the partial derivative of T with respect to theta2
    T1_wrt_theta2 = BuildT1(Partial_fw_wrt_theta2,n1-t1,n2-t2);
    T2_wrt_theta2 = BuildT1(Partial_gw_wrt_theta2,m1-t1,m2-t2);
    
    Partial_T_wrt_theta2 = [T1_wrt_theta2 alpha(ite)*T2_wrt_theta2];
    
    % Calculate the column c_{k} of DTQ that is moved to the right hand side
    ck = T*e;
    
    % Calculate the derivatives of c_{k} with respect to \alpha and \theta
    Partial_ck_wrt_alpha        = Partial_T_wrt_alpha*e;
    Partial_ck_wrt_theta1       = Partial_T_wrt_theta1*e;
    Partial_ck_wrt_theta2       = Partial_T_wrt_theta2*e;
    
    % Create the vector of structured perturbations zf and zg applied
    % to F and G.
    z_fx      = zk(1:num_coeff_f);
    z_gx      = zk(num_coeff_f + 1 :end);
    
    % Calculate the entries of z_fw and z_gw
    z_fx_mat = GetAsMatrix(z_fx,m1,m2);
    z_gx_mat = GetAsMatrix(z_gx,n1,n2);
    
    % Get z_fw_mat, by multiplying by thetas
    z_fw_mat = GetWithThetas(z_fx_mat,th1(ite),th2(ite));
    z_gw_mat = GetWithThetas(z_gx_mat,th1(ite),th2(ite));
    
    % Calculate the derivatives of z_fw and z_gw with repect to alpha.
    Partial_zfw_wrt_alpha    = zeros(m1+1,m2+1);
    Partial_zgw_wrt_alpha    = z_gw_mat;
    
    % Calculate the derivative of z_fw with respect to theta1.
    temp_mat = diag((0:1:m1)./th1(ite));
    Partial_zfw_wrt_theta1 = temp_mat * z_fw_mat;
    
    % Calculate the derivative of z_fw with respect to theta2
    temp_mat = diag((0:1:m2)./th2(ite));
    Partial_zfw_wrt_theta2 = z_fw_mat * temp_mat;
    
    % Calculate the derivative of z_gw with respect ot theta1
    temp_mat = diag((0:1:n1)./th2(ite));
    Partial_zgw_wrt_theta1 = temp_mat * z_gw_mat;
    
    % Calculate the deriviate of z_gw with respect to theta2
    temp_mat = diag((0:1:n2)./th2(ite));
    Partial_zgw_wrt_theta2 = z_gw_mat * temp_mat;
    
    % Build the coefficient Matrix N = [T(z1) T(z2)], of structured perturbations, with
    % same structure as DTQ.
    N1 = BuildT1(z_fw_mat,n1-t1,n2-t2);
    N2 = BuildT1(z_gw_mat,m1-t1,m2-t2);
    N = [N1 alpha(ite).*N2];
    
    % Build the coefficient matrix N with respect to alpha
    T1_wrt_alpha = BuildT1(Partial_zfw_wrt_alpha, n1-t1,n2-t2);
    T2_wrt_alpha = BuildT1(Partial_zgw_wrt_alpha, m1-t1,m2-t2);
    Partial_N_wrt_alpha = [T1_wrt_alpha T2_wrt_alpha];
    
    
    % Calculate the derivatives of DNQ with respect to theta
    T1_wrt_theta1 = BuildT1(Partial_zfw_wrt_theta1,n1-t1,n2-t2);
    T2_wrt_theta1 = BuildT1(Partial_zgw_wrt_theta1,m1-t1,m2-t2);
    Partial_N_wrt_theta1 = [T1_wrt_theta1 alpha(ite).*T2_wrt_theta1];
    
    % Calculate the derivatives of DNQ with respect to theta
    T1_wrt_theta2 = BuildT1(Partial_zfw_wrt_theta2,n1-t1,n2-t2);
    T2_wrt_theta2 = BuildT1(Partial_zgw_wrt_theta2,m1-t1,m2-t2);
    Partial_N_wrt_theta2 = [T1_wrt_theta2 alpha(ite).*T2_wrt_theta2];
    
    % Calculate the column of DNQ that is moved to the right hand side, which
    % has the same structure as c_{k} the column of S_{k} moved to the RHS
    h = N*e;
    
    % Calculate the derivative of h with respect to alpha
    h_alpha = Partial_N_wrt_alpha*e;
    
    % Calculate the derivative of h with respect to theta1
    h_theta1 = Partial_N_wrt_theta1*e;
    
    % Calculate the derivative of h with respect to theta2
    h_theta2 = Partial_N_wrt_theta2*e;
    
    % Build the matrix (T+N)
    
    T1 = BuildT1(fww_matrix + z_fw_mat,n1-t1,n2-t2);
    T2 = BuildT1(gww_matrix + z_gw_mat,m1-t1,m2-t2);
    TN = [T1 alpha(ite).*T2];
    
    % Calculate the paritial derivative of (T+N) with respect to
    % alpha
    TN1_wrt_alpha = BuildT1(Partial_fw_wrt_alpha + Partial_zfw_wrt_alpha, n1-t1,n2-t2);
    TN2_wrt_alpha = BuildT1(Partial_alpha_gw_wrt_alpha + Partial_zgw_wrt_alpha, m1-t1,m2-t2);
    TN_alpha = [TN1_wrt_alpha TN2_wrt_alpha];
    
    
    % Calculate the paritial derivative of (T+N) with respect to
    % theta1
    TN1_wrt_theta1 = BuildT1(Partial_fw_wrt_theta1 + Partial_zfw_wrt_theta1,n1-t1,n2-t2);
    TN2_wrt_theta1 = BuildT1(Partial_gw_wrt_theta1 + Partial_zgw_wrt_theta1,m1-t1,m2-t2);
    TN_theta1 = [TN1_wrt_theta1 alpha(ite).*TN2_wrt_theta1];
    
    % Calculate the paritial derivative of (T+N) with respect to
    % theta2
    TN1_wrt_theta2 = BuildT1(Partial_fw_wrt_theta2 + Partial_zfw_wrt_theta2,n1-t1,n2-t2);
    TN2_wrt_theta2 = BuildT1(Partial_gw_wrt_theta2 + Partial_zgw_wrt_theta2,m1-t1,m2-t2);
    TN_theta2 = [TN1_wrt_theta2 alpha(ite).*TN2_wrt_theta2];
    
    % Update xk
    xk = SolveAx_b(TN*M,ck+h);
    
    % Calculate the matrix DY where Y is the Matrix such that E_{k}x = Y_{k}z.
    Y = BuildY(m1,m2,n1,n2,t1,t2,opt_col,xk,alpha(ite),th1(ite),th2(ite));
    
    % Calculate the matrix P where ck = P * [f,g]
    P = BuildP(m1,m2,n1,n2,alpha(ite),th1(ite),th2(ite),opt_col,t1,t2,num_cols_T1);
    
    % Test P
    %ck;
    %ck2 = P*[fxy_vec;gxy_vec];
    %ck - ck2;
    
    
    %%
    
    % Get residual as a vector
    rk = (ck+h) - TN*M*xk ;
    
    % Create the matrix C. This is made up of five submatrices, HZ, Hx,
    % H_alpha and H_theta1 and H_theta2.
    
    Hz          = Y-P;
    
    Hx          = TN*M;
    
    H_alpha     = TN_alpha*M*xk - (Partial_ck_wrt_alpha + h_alpha);
    
    H_theta1    = TN_theta1*M*xk - (Partial_ck_wrt_theta1 + h_theta1);
    
    H_theta2    = TN_theta2*M*xk - (Partial_ck_wrt_theta2 + h_theta2);
    
    C = [Hz,Hx,H_alpha,H_theta1, H_theta2];  % the matrix C
    
    % Calculate the new right hand vector
    ek = ck + h;
    
    % update gnew - used in LSE Problem.
    res_vec = rk;
    
    % Calculate the normalised residual of the solution.
    condition(ite) = norm(rk) / norm(ek);
    
    % Update fnew - used in LSE Problem.
    f = -(yy-start_point);
    
end

%% Plot Graphs
switch PLOT_GRAPHS
    case 'y'
        figure('name','Residuals in SNTLN')
        hold on
        title('Residuals in SNTLN')
        xlabel('Iterations')
        ylabel('log_{10} Residuals')
        plot((1:1:ite),log10(condition),'-s','DisplayName','Residual')
        legend(gca,'show');
        hold off
        
        figure('name','Theta Variation over Newton Raphson Iterations')
        hold on
        title('Variation of \theta over Newton Raphson Iteration')
        plot((1:1:ite),log10(th1),'-s','DisplayName','\theta_{1}')
        plot((1:1:ite),log10(th2),'-s','DisplayName','\theta_{2}')
        legend(gca,'show');
        hold off
        
        figure('name','Alpha variation over Newton Raphson Iterations')
        hold on
        title('Alpha variation over Newton Raphson Iterations')
        plot((1:1:ite),log10(alpha),'-s','DisplayName','\alpha')
        legend(gca,'show');
        hold off
    case 'n'
    otherwise
        error('plot_graphs is either y or n')
end


%%

% Once iterations are complete, assign fx output, gx output, solution X
% output, alpha output and theta output.

% get the vector zk
zPert_f_vec = zk(1:num_coeff_f);
zPert_f_mat = GetAsMatrix(zPert_f_vec,m1,m2);

zPert_g_vec = zk(num_coeff_f+1:end);
zPert_g_mat = GetAsMatrix(zPert_g_vec,n1,n2);

% Set outputs of low rank approximation

o_fxy = fxy_matrix + zPert_f_mat;

o_gxy = gxy_matrix + zPert_g_mat;

o_X  = xk;

o_alpha = alpha(ite);

o_th1 = th1(ite);

o_th2 = th2(ite);


if ite == MAX_ITERATIONS_SNTLN
    
    
    if condition(ite) < condition(1)
        val = 'Keep Final'
    else
        val = 'Revert'
    end
    switch val
        case 'Revert'
            fprintf('SNTLN Failed to converge, default to input values\n')
            % SNTLN Failed so revert to previous values
            o_fxy = fxy_matrix;
            o_gxy = gxy_matrix;
            o_alpha = i_alpha;
            o_th1 = i_th1;
            o_th2 = i_th2;
            o_X = initial_xls;
            return;
        case 'Keep Final'
            fprintf('SNTLN Failed to converge, keep termination values\n')
            o_fxy = fxy_matrix + zPert_f_mat;
            o_gxy = gxy_matrix + zPert_g_mat;
            o_X  = xk;
            o_alpha = alpha(ite);
            o_th1 = th1(ite);
            o_th2 = th2(ite);
            
            return;
    end
    
end


% Print the number of iterations
LineBreakLarge();
fprintf([mfilename sprintf('SNTLN converged within %i iterations \n', ite)]);
LineBreakLarge();

end

















