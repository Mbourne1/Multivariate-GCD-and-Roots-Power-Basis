function [ fxy_lr, gxy_lr, alpha_lr,th1_lr,th2_lr,X_lr] = ...
    SNTLN_Total(fxy_matrix, gxy_matrix, m, n, i_alpha, i_th1, i_th2,k,idx_col)
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
% fxy_lr : Coefficients of f(x,y) on output, in standard bernstein basis,
% including added structured perturbations.
%
% gxy_lr : Coefficients of fx on output, in standard bernstein basis,
% including added structured perturbations.
%
% alpha_lr : Optimal value of \alpha
%
% th1_lr : Optimal value of \theta_{1}
%
% th2_lr : Optimal value of \theta_{2}
%
% X_lr :
%

% Global Inputs
global SETTINGS

%%
% Pad the coefficients of fxy and gxy
% this is equivalent to degree elevating so that f is of degree (m,m), and
% g is of degree (n,n)
fxy_matrix_padd = zeros(m+1,m+1);
gxy_matrix_padd = zeros(n+1,n+1);


[r,c] = size(fxy_matrix);
fxy_matrix_padd(1:r,1:c) = fxy_matrix;

[r,c] = size(gxy_matrix);
gxy_matrix_padd(1:r,1:c) = gxy_matrix;

fxy_matrix = fxy_matrix_padd;
gxy_matrix = gxy_matrix_padd;

%%

% Set the initial iterations number
ite = 1;

% Set initial values of alpha and theta
th1(ite) = i_th1;
th2(ite) = i_th2;
alpha(ite) = i_alpha;

% Get the number of coefficients in the polynomial f(x,y)
nCoeff_fxy = (m+1) * (m+1);
nNonZeros_fxy = nchoosek(m+2,2);
nZeros_fxy = nchoosek(m+1,2);

% Get the number of coefficients in the polynomial g(x,y)
nCoeff_gxy = (n+1) * (n+1);
nNonZeros_gxy = nchoosek(n+2,2);
nZeros_gxy = nchoosek(n+1,2);

% Get the number of coefficients in both f(x,y) and g(x,y)
nCoeff_fg = nCoeff_fxy + nCoeff_gxy;
nNonZeros_fg = nNonZeros_fxy + nNonZeros_gxy;
nZeros_fg = nZeros_fxy + nZeros_gxy;

% Get the number of coefficients in v(x,y)
nCoeff_vxy = (n-k+1) * (n-k+1);
nNonZeros_vxy = nchoosek(n-k+2,2);
nZeros_vxy = nchoosek(n-k+1,2);

% Get the number of coefficients in u(x,y)
nCoeff_uxy = (m-k+1) * (m-k+1);
nNonZeros_uxy = nchoosek(m-k+2,2);
nZeros_uxy = nchoosek(m-k+1,2);

% Get the number of coefficients in the unknown vector x, where A_{t}x =
% c_{t}.
nCoeff_x = nCoeff_uxy + nCoeff_vxy - 1;
nNonZeros_x = nNonZeros_uxy + nNonZeros_vxy - 1;

% Create the identity matrix I, the matrix M formed from I by removing the
% column equivalent to the optimal column for removal from the Sylvester
% subresultant matrix, so that S_{t}(f,g)*M = A_{t}, where A_{t} is the
% Sylvester subresultant matrix with the column removed.

% Get the number of columns in C_{t}(f), the first partition of the Sylvester
% Matrix S_{t}(f,g)
nCols_T1 = nNonZeros_vxy;

% Get the number of columns in C_{t}(g), the second partition of the
% Sylvester matrix S_{t}(f,g)
nCols_T2 = nNonZeros_uxy;

% Get the total number of columns in the Sylvester matrix S_{k}(f,g)
nCols_Sk = nCols_T1 + nCols_T2;

% Get the number of rows in the Sylvester subresultant matrix S_{k}(f,g)
nRows_Sk = nchoosek(m+n-k+2,2);


% Create the identity matrix
I = eye(nCols_Sk, nCols_Sk);

% Create the matrix M, such that S(f,g)*M gives A_{t}, the Sylvester Matrix
% with the optimal column removed.
M = I;
M(:,idx_col) = [];

% Let e be the column removed from the identity matrix, such that
% S_{t}(f,g) * e gives the column c_{t}, where c_{t} is the optimal column
% removed from the Sylvester subresultant.
e = I(:,idx_col);

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
T1_fxy = BuildT1_Total(fww_matrix,m,n-k);
T2_gxy = BuildT1_Total(gww_matrix,n,m-k);
T_fg = [T1_fxy alpha(ite).*T2_gxy];

%
% Calculate the partial derivatives of f(w,w) and g(w,w) with respect to \alpha
Partial_fw_wrt_alpha            = zeros(m+1,m+1);
Partial_alpha_gw_wrt_alpha      = gxy_matrix;

%
% Calculate the partial derivatives of f(w,w) with respect to \theta_1
theta_mat = diag((0:1:m) ./ th1(ite));
Partial_fw_wrt_theta1    = theta_mat * fww_matrix;

%
% Calculate the partial derivative of g(w,w) with respect to theta_1
theta_mat = diag((0:1:n) ./ th1(ite));
Partial_gw_wrt_theta1    = theta_mat * gww_matrix;

%
% Calculate the partial derivative of f(w,w) with respect to theta_2
theta_mat = diag((0:1:m) ./ th2(ite));
Partial_fw_wrt_theta2 = fww_matrix * theta_mat;

%
% Calculate the partial deriviates of g(w,w) with respect to theta_2
theta_mat = diag((0:1:n)./ th2(ite));
Partial_gw_wrt_theta2 = gww_matrix * theta_mat;

%%
% Build the derivative of T(f,g) with respect to alpha
T1_f_wrt_alpha = BuildT1_Total(Partial_fw_wrt_alpha,m,n-k);
T2_g_wrt_alpha = BuildT1_Total(Partial_alpha_gw_wrt_alpha,n,m-k);
T_fg_wrt_alpha = [T1_f_wrt_alpha T2_g_wrt_alpha];

%%
% Calculate the derivative of T(f,g) with respect to theta_{1}
T1_f_wrt_theta1 = BuildT1_Total(Partial_fw_wrt_theta1,m,n-k);
T2_g_wrt_theta1 = BuildT1_Total(Partial_gw_wrt_theta1,n,m-k);
Partial_T_fg_wrt_theta1 = [T1_f_wrt_theta1 alpha(ite)* T2_g_wrt_theta1];

%%
% Calcualte the derivative of T(f,g) with respect to theta_2
T1_f_wrt_theta2 = BuildT1_Total(Partial_fw_wrt_theta2,m,n-k);
T2_g_wrt_theta2 = BuildT1_Total(Partial_gw_wrt_theta2,n,m-k);
Partial_T_fg_wrt_theta2 = [T1_f_wrt_theta2 alpha(ite)*T2_g_wrt_theta2];

%%
% Initialise the vector z of structured perturbations
% if we are working with strictly the roots problem, the number of entries
% in z can be reduced.


zk = zeros(nNonZeros_fg , 1);

%%
% Initilaise the derivative of N wrt alpha.
Partial_N_wrt_alpha   = zeros(nRows_Sk,nCols_T1 + nCols_T2);

% Initilaise the derivative of N wrt theta_1.
Partial_N_wrt_theta1   = zeros(nRows_Sk,nCols_T1 + nCols_T2);

% Initialise the derivative of N wrt theta 2
Partial_N_wrt_theta2   = zeros(nRows_Sk,nCols_T1 + nCols_T2);

%%
% Initialise the derivative of h
% Calculate the derivatives wrt alpha and theta of the column of DNQ
% that is moved to the right hand side.
Partial_h_wrt_alpha     = Partial_N_wrt_alpha*e;
Partial_h_wrt_theta1    = Partial_N_wrt_theta1*e;
Partial_h_wrt_theta2    = Partial_N_wrt_theta2*e;


% Get the matrix A_{k}(f,g), which is the subresultant matrix S(f,g) with
% an opitmal column removed
Ak = T_fg;
ck = T_fg(:,idx_col);
Ak(:,idx_col) = [];

%% Build the matrix P
P = BuildP_TotalDegree_SNTLN(m,n,k,alpha(ite),th1(ite),th2(ite),idx_col);
% Test P
% Get the coefficients of f(x,y) in matrix form
%%fxy_vec = GetAsVector(fxy_matrix);
%%gxy_vec = GetAsVector(gxy_matrix);

%%ck;
%%ck2 = P*[fxy_vec;gxy_vec];
%%ck - ck2;

%%
% Calculate the derivatives of ck wrt alpha and theta.
Partial_ck_wrt_alpha        = T_fg_wrt_alpha*e;
Partial_ck_wrt_theta1       = Partial_T_fg_wrt_theta1*e;
Partial_ck_wrt_theta2       = Partial_T_fg_wrt_theta2*e;


%%
% Perform QR decomposition of Ak to obtain the solution x
x_ls = SolveAx_b(Ak,ck);

% store this version of x_ls
initial_xls = x_ls;
first_part = x_ls(1:(idx_col-1));
second_part = x_ls(idx_col:end);
x = [first_part ; 0 ; second_part];

% Build Matrix Y, where Y(v,u)*[f;g] = S(f,g)*[u;v]
Y = BuildY_TotalDegree_SNTLN(x,m,n,k,alpha(ite),th1(ite),th2(ite));
% Test Y
%test1 = Y * [fxy_vec;gxy_vec];
%test2 = T * x;
%test1-test2;

% Calculate the initial residual r = ck - (Ak*x)
res_vec = ck - (T_fg*M*x_ls);

% Get the matrix p, which will store all the perturbations returned from LSE file
nEntries = nNonZeros_fxy...
    + nNonZeros_gxy ...
    + nNonZeros_x ...
    + 3;


% Set the initial value of vector p to be zero
f = zeros(nNonZeros_fxy...
    + nNonZeros_gxy ...
    + nNonZeros_x ...
    + 3,1);

%
% Set the intial value of E to the identity matrix
E = eye(nEntries);

%
% Create the matrix (T+N), initially N is empty so this is the same as T.
T1_fxy = BuildT1_Total(fww_matrix,m,n-k);
T2_gxy = BuildT1_Total(gww_matrix,n,m-k);
TN = [T1_fxy alpha(ite).*T2_gxy];

%
% Create The matrix (T+N) with respect to alpha
T1_f_wrt_alpha = BuildT1_Total(Partial_fw_wrt_alpha,m,n-k);
T2_g_wrt_alpha = BuildT1_Total(Partial_alpha_gw_wrt_alpha,n,m-k);
TN_wrt_alpha = [T1_f_wrt_alpha T2_g_wrt_alpha];

%
% Create The matrix (T+N) with respect to theta1
T1_f_wrt_theta1 = BuildT1_Total(Partial_fw_wrt_theta1,m,n-k);
T2_g_wrt_theta1 = BuildT1_Total(Partial_gw_wrt_theta1,n,m-k);
TN_wrt_theta1 = [T1_f_wrt_theta1 alpha(ite) * T2_g_wrt_theta1];

%
% Create The matrix (T+N) with respect to theta2
T1_f_wrt_theta2 = BuildT1_Total(Partial_fw_wrt_theta2,m,n-k);
T2_g_wrt_theta2 = BuildT1_Total(Partial_gw_wrt_theta2,n,m-k);
TN_wrt_theta2 = [T1_f_wrt_theta2 alpha(ite) * T2_g_wrt_theta2];

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

yy =  start_point;

% Set the termination criterion to a large value. It will be
% over written later.
condition(ite) = norm(res_vec) ./ norm(ck);

xk = x_ls;


while condition(ite) >(SETTINGS.MAX_ERROR_SNTLN) &&  ite < SETTINGS.MAX_ITERATIONS_SNTLN
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
    
    % Get the coefficients corresponding to f and g
    delta_zk        = y(1:nNonZeros_fxy + nNonZeros_gxy ,1);
    
    % Remove the zk coefficients from the list of coefficients
    y(1:nNonZeros_fxy + nNonZeros_gxy) = [];
    
    % Get the coefficients corresponding to x
    delta_xk        = y(1:nNonZeros_x,1);
    
    % Remove them from the list of coefficients
    y(1:nNonZeros_x) = [];
    
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
    T1_fxy = BuildT1_Total(fww_matrix,m,n-k);
    T2_gxy = BuildT1_Total(gww_matrix,n,m-k);
    T_fg = [T1_fxy alpha(ite).*T2_gxy];
    
    % Calculate the partial derivatives of fw and gw with respect to alpha
    Partial_fw_wrt_alpha            = zeros(m+1,m+1);
    Partial_alpha_gw_wrt_alpha      = gww_matrix;
    
    %%
    % Calculate the partial derivatives of fw and gw with respect to theta1
    % divide the rows by theta1 and multiply by the old power
    
    % Get the partial derivative of f with respect to theta 1
    temp_mat = diag((0:1:m)./th1(ite));
    Partial_fw_wrt_theta1 = temp_mat * fww_matrix;
    
    % Get the partial derivative of g with respect to theta1
    temp_mat = diag((0:1:n)./th1(ite));
    Partial_gw_wrt_theta1 = temp_mat * gww_matrix;
    
    % Get the partial derivative of f with respect to theta2
    temp_mat = diag((0:1:m)./th2(ite));
    Partial_fw_wrt_theta2 =  fww_matrix * temp_mat;
    
    % Get the partial derivative of g with respect too theta2
    temp_mat = diag((0:1:n)./th2(ite));
    Partial_gw_wrt_theta2 =  gww_matrix * temp_mat;
    
    % Calculate the Partial derivative of T with respect to alpha.
    T1_f_wrt_alpha = BuildT1_Total(Partial_fw_wrt_alpha, m, n-k);
    T2_g_wrt_alpha = BuildT1_Total(Partial_alpha_gw_wrt_alpha, n, m-k);
    Partial_T_wrt_alpha = [T1_f_wrt_alpha T2_g_wrt_alpha];
    
    % Calculate the partial derivative of T with respect to theta1
    T1_f_wrt_theta1 = BuildT1_Total(Partial_fw_wrt_theta1,m,n-k);
    T2_g_wrt_theta1 = BuildT1_Total(Partial_gw_wrt_theta1,n,m-k);
    Partial_T_fg_wrt_theta1 = [T1_f_wrt_theta1 alpha(ite)*T2_g_wrt_theta1];
    
    % Calculate the partial derivative of T with respect to theta2
    T1_f_wrt_theta2 = BuildT1_Total(Partial_fw_wrt_theta2,m,n-k);
    T2_g_wrt_theta2 = BuildT1_Total(Partial_gw_wrt_theta2,n,m-k);
    
    Partial_T_fg_wrt_theta2 = [T1_f_wrt_theta2 alpha(ite)*T2_g_wrt_theta2];
    
    % Calculate the column c_{k} of DTQ that is moved to the right hand side
    ck = T_fg*e;
    
    % Calculate the derivatives of c_{k} with respect to \alpha and \theta
    Partial_ck_wrt_alpha        = Partial_T_wrt_alpha*e;
    Partial_ck_wrt_theta1       = Partial_T_fg_wrt_theta1*e;
    Partial_ck_wrt_theta2       = Partial_T_fg_wrt_theta2*e;
    
    % Create the vector of structured perturbations zf and zg applied
    % to F and G.
    z_fx      = zk(1:nNonZeros_fxy);
    z_gx      = zk(nNonZeros_fxy + 1 :end);
    
    % Calculate the entries of z_fw and z_gw
    z_fx_mat = GetAsMatrix([z_fx;zeros(nZeros_fxy,1)],m,m);
    z_gx_mat = GetAsMatrix([z_gx;zeros(nZeros_gxy,1)],n,n);
    
    % Get z_fw_mat, by multiplying by thetas
    z_fw_mat = GetWithThetas(z_fx_mat,th1(ite),th2(ite));
    z_gw_mat = GetWithThetas(z_gx_mat,th1(ite),th2(ite));
    
    % Calculate the derivatives of z_fw and z_gw with repect to alpha.
    Partial_zfw_wrt_alpha    = zeros(m+1,m+1);
    Partial_zgw_wrt_alpha    = z_gw_mat;
    
    % Calculate the derivative of z_fw with respect to theta1.
    temp_mat = diag((0:1:m)./th1(ite));
    Partial_zfw_wrt_theta1 = temp_mat * z_fw_mat;
    
    % Calculate the derivative of z_fw with respect to theta2
    temp_mat = diag((0:1:m)./th2(ite));
    Partial_zfw_wrt_theta2 = z_fw_mat * temp_mat;
    
    % Calculate the derivative of z_gw with respect ot theta1
    temp_mat = diag((0:1:n)./th2(ite));
    Partial_zgw_wrt_theta1 = temp_mat * z_gw_mat;
    
    % Calculate the deriviate of z_gw with respect to theta2
    temp_mat = diag((0:1:n)./th2(ite));
    Partial_zgw_wrt_theta2 = z_gw_mat * temp_mat;
    
    % Build the coefficient Matrix N = [T(z1) T(z2)], of structured perturbations, with
    % same structure as DTQ.
    N1 = BuildT1_Total(z_fw_mat,m,n-k);
    N2 = BuildT1_Total(z_gw_mat,n,m-k);
    N = [N1 alpha(ite).*N2];
    
    % Build the coefficient matrix N with respect to alpha
    T1_f_wrt_alpha = BuildT1_Total(Partial_zfw_wrt_alpha, m, n-k);
    T2_g_wrt_alpha = BuildT1_Total(Partial_zgw_wrt_alpha, n, m-k);
    Partial_N_wrt_alpha = [T1_f_wrt_alpha T2_g_wrt_alpha];
    
    
    % Calculate the derivatives of DNQ with respect to theta
    T1_f_wrt_theta1 = BuildT1_Total(Partial_zfw_wrt_theta1, m, n-k);
    T2_g_wrt_theta1 = BuildT1_Total(Partial_zgw_wrt_theta1, n, m-k);
    Partial_N_wrt_theta1 = [T1_f_wrt_theta1 alpha(ite).*T2_g_wrt_theta1];
    
    % Calculate the derivatives of DNQ with respect to theta
    T1_f_wrt_theta2 = BuildT1_Total(Partial_zfw_wrt_theta2, m, n-k);
    T2_g_wrt_theta2 = BuildT1_Total(Partial_zgw_wrt_theta2, n, m-k);
    Partial_N_wrt_theta2 = [T1_f_wrt_theta2 alpha(ite).*T2_g_wrt_theta2];
    
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
    
    T1_fxy = BuildT1_Total(fww_matrix + z_fw_mat, m, n-k);
    T2_gxy = BuildT1_Total(gww_matrix + z_gw_mat, n, m-k);
    TN = [T1_fxy alpha(ite).*T2_gxy];
    
    % Calculate the paritial derivative of (T+N) with respect to
    % alpha
    TN1_wrt_alpha = BuildT1_Total(Partial_fw_wrt_alpha + Partial_zfw_wrt_alpha, m, n-k);
    TN2_wrt_alpha = BuildT1_Total(Partial_alpha_gw_wrt_alpha + Partial_zgw_wrt_alpha, n, m-k);
    TN_alpha = [TN1_wrt_alpha TN2_wrt_alpha];
    
    
    % Calculate the paritial derivative of (T+N) with respect to
    % theta1
    TN1_wrt_theta1 = BuildT1_Total(Partial_fw_wrt_theta1 + Partial_zfw_wrt_theta1, m, n-k);
    TN2_wrt_theta1 = BuildT1_Total(Partial_gw_wrt_theta1 + Partial_zgw_wrt_theta1, n, m-k);
    TN_theta1 = [TN1_wrt_theta1 alpha(ite).*TN2_wrt_theta1];
    
    % Calculate the paritial derivative of (T+N) with respect to
    % theta2
    TN1_wrt_theta2 = BuildT1_Total(Partial_fw_wrt_theta2 + Partial_zfw_wrt_theta2, m, n-k);
    TN2_wrt_theta2 = BuildT1_Total(Partial_gw_wrt_theta2 + Partial_zgw_wrt_theta2, n, m-k);
    TN_theta2 = [TN1_wrt_theta2 alpha(ite).*TN2_wrt_theta2];
    
    
    % Calculate the matrix DY where Y is the Matrix such that E_{k}x = Y_{k}z.
    first_part = xk(1:(idx_col-1));
    second_part = xk(idx_col:end);
    x = [first_part ; 0 ; second_part];
    
    Y = BuildY_TotalDegree_SNTLN(x,m,n,k,alpha(ite),th1(ite),th2(ite));
    
    % Calculate the matrix P where ck = P * [f,g]
    P = BuildP_TotalDegree_SNTLN(m,n,k,alpha(ite),th1(ite),th2(ite),idx_col);
    
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

PlotGraphs_SNTLN()


%

% Once iterations are complete, assign fx output, gx output, solution X
% output, alpha output and theta output.

% get the vector zk
zPert_f_vec = zk(1:nNonZeros_fxy);
zPert_f_mat = GetAsMatrix([zPert_f_vec; zeros(nZeros_fxy,1)],m,m);

zPert_g_vec = zk(nNonZeros_fxy+1:end);
zPert_g_mat = GetAsMatrix([zPert_g_vec; zeros(nZeros_gxy,1)],n,n);

% Set outputs of low rank approximation

fprintf([mfilename ' : ' sprintf('SNTLN Failed to converge, keep termination values\n')])
fxy_lr = fxy_matrix + zPert_f_mat;
gxy_lr = gxy_matrix + zPert_g_mat;
X_lr  = xk;
alpha_lr = alpha(ite);
th1_lr = th1(ite);
th2_lr = th2(ite);


% Print the number of iterations
LineBreakLarge();
fprintf([mfilename ' : ' sprintf('SNTLN converged within %i iterations \n', ite)]);
LineBreakLarge()

end

















