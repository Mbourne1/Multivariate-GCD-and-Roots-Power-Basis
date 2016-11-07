function [ fxy_lr, gxy_lr, alpha_lr,th1_lr,th2_lr,X_lr] = ...
    SNTLN_Both( fxy_matrix,gxy_matrix, i_alpha, i_th1, i_th2,m,n,k,k1,k2,idx_col)
% Obtain the low rank approximation of the Sylvester matrix D*T_{t}(f,g)*Q =
% S_{t}(f,g)
%
% SNTLN( fxy_matrix,gxy_matrix, i_alpha, i_th1, i_th2,t1,t2,opt_col)
%
% % Inputs:
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
% % Outputs:
%
% o_fxy : Coefficients of f(x,y) on output, in standard bernstein basis,
% including added structured perturbations.
%
% o_gxy : Coefficients of g(x,y) on output, in standard bernstein basis,
% including added structured perturbations.
%
% o_alpha : Optimal value of \alpha
%
% o_th1 : Optimal value of \theta_{1}
%
% o_th2 : Optimal value of \theta_{2}
%
% o_X : 
%

% Global Inputs

global SETTINGS


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
nCoeff_fxy =  (m1+1) * (m2+1);
nNonZeros_fxy = GetNumNonZeros(m1,m2,m);
nZeros_fxy = nCoeff_fxy - nNonZeros_fxy;

% Get the number of coefficients in the polynomial g(x,y)
nCoeff_gxy = (n1+1) * (n2+1);
nNonZeros_gxy = GetNumNonZeros(n1,n2,n);
nZeros_gxy = nCoeff_gxy - nNonZeros_gxy;

% 
nCoeff_fv = GetNumNonZeros(m1+n1-k1,m2+n2-k2,m+n-k);

% Get the number of coefficients in v(x,y)
nCoeff_vxy = (n1-k1+1) * (n2-k2+1);
nNonZeros_vxy = GetNumNonZeros(n1-k1,n2-k2,n-k);
nZeros_vxy = nCoeff_vxy - nNonZeros_vxy;

% Get the number of coefficients in u(x,y)
nCoeff_uxy = (m1-k1+1) * (m2-k2+1);
nNonZeros_uxy = GetNumNonZeros(m1-k1,m2-k2,m-k);
nZeros_uxy = nCoeff_uxy - nNonZeros_uxy;

% Get the number of coefficients in the product f*v(x,y)

% Get the number of coefficients in the unknown vector x, where A_{t}x =
% c_{t}.
nCoeff_x = nNonZeros_uxy + nNonZeros_vxy - 1;

% Create the identity matrix I, the matrix M formed from I by removing the
% column equivalent to the optimal column for removal from the Sylvester
% subresultant matrix, so that S_{t}(f,g)*M = A_{t}, where A_{t} is the
% Sylvester subresultant matrix with the column removed.

% Get the number of columns in C_{t}(f), the first partition of the Sylvester
% Matrix S_{t}(f,g)
nCols_Tf = nNonZeros_vxy;


% Get the number of columns in C_{t}(g), the second partition of the
% Sylvester matrix S_{t}(f,g)
nCols_Tg = nNonZeros_uxy;

% Get the total number of columns in the Sylvester matrix S_{t}(f,g)
nCols_Sylv = nCols_Tf + nCols_Tg;

% Create the identity matrix
I = eye(nCols_Sylv, nCols_Sylv);

% Create the matrix M, such that S(f,g)*M gives A_{t}, the Sylvester Matrix
% with the optimal column removed.
M = I;
M(:,idx_col) = [];

% Let e be the column removed from the identity matrix, such that
% S_{t}(f,g) * e gives the column c_{t}, where c_{t} is the optimal column
% removed from the Sylvester subresultant.
e = I(:,idx_col);

% % Preprocessing

% Obtain polynomials in Modified Bernstein Basis, using initial values of
% alpha and theta.

% Multiply the rows of fxy_matrix by theta1, and multiply the cols of
% fxy_matrix by theta2
fww_matrix = GetWithThetas(fxy_matrix,th1(ite),th2(ite));

% Multiply the rows of gxy_matrix by theta1, and multiply the cols of
% gxy_matrix by theta2
gww_matrix = GetWithThetas(gxy_matrix,th1(ite),th2(ite));

% Form the Coefficient Matrix T = [C(f(w1,w2))|alpha * C(g(w1,w2))] such that T*x = [col]
Tf = BuildT1_Both(fww_matrix,m,n-k,n1-k1,n2-k2);
Tg = BuildT1_Both(gww_matrix,n,m-k,m1-k1,m2-k2);
Sk_fg = [Tf alpha(ite).*Tg];

%
% Calculate the partial derivatives of f(w,w) and g(w,w) with respect to \alpha
Partial_fw_wrt_alpha_mat            = zeros(m1+1,m2+1);
Partial_alpha_gw_wrt_alpha_mat      = gxy_matrix;

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

% %
% Build the derivative of T(f,g) with respect to alpha
Tf_wrt_alpha = BuildT1_Both(Partial_fw_wrt_alpha_mat,m,n-k,n1-k1,n2-k2);
Tg_wrt_alpha = BuildT1_Both(Partial_alpha_gw_wrt_alpha_mat,n,m-k,m1-k1,m2-k2);
T_alpha = [Tf_wrt_alpha Tg_wrt_alpha];

% %
% Calculate the derivative of T(f,g) with respect to theta_{1}
Tf_wrt_theta1 = BuildT1_Both(Partial_fw_wrt_theta1,m,n-k,n1-k1,n2-k2);
Tg_wrt_theta1 = BuildT1_Both(Partial_gw_wrt_theta1,n,m-k,m1-k1,m2-k2);
Partial_T_wrt_theta1 = [Tf_wrt_theta1 alpha(ite)* Tg_wrt_theta1];

% %
% Calcualte the derivative of T(f,g) with respect to theta_2
Tf_wrt_theta2 = BuildT1_Both(Partial_fw_wrt_theta2,m,n-k,n1-k1,n2-k2);
Tg_wrt_theta2 = BuildT1_Both(Partial_gw_wrt_theta2,n,m-k,m1-k1,m2-k2);
Partial_T_wrt_theta2 = [Tf_wrt_theta2 alpha(ite)*Tg_wrt_theta2];

% %
% Initialise the vector z of structured perturbations
% if we are working with strictly the roots problem, the number of entries
% in z can be reduced.

nRows_Sylv_mat = nCoeff_fv;
nCols_Tf = nNonZeros_vxy;
nCols_Tg = nNonZeros_uxy;

zk = zeros(nNonZeros_fxy + nNonZeros_gxy , 1);

%
% Initilaise the derivative of N wrt alpha.
Partial_N_wrt_alpha   = zeros(nRows_Sylv_mat,nCols_Tf + nCols_Tg);

% Initilaise the derivative of N wrt theta_1.
Partial_N_wrt_theta1   = zeros(nRows_Sylv_mat,nCols_Tf + nCols_Tg);

% Initialise the derivative of N wrt theta 2
Partial_N_wrt_theta2   = zeros(nRows_Sylv_mat,nCols_Tf + nCols_Tg);

%
% Initialise the derivative of h
% Calculate the derivatives wrt alpha and theta of the column of DNQ
% that is moved to the right hand side.
Partial_h_wrt_alpha     = Partial_N_wrt_alpha*e;
Partial_h_wrt_theta1    = Partial_N_wrt_theta1*e;
Partial_h_wrt_theta2    = Partial_N_wrt_theta2*e;


% Get the matrix A_{k}(f,g), which is the subresultant matrix S(f,g) with
% an opitmal column removed
Ak_fg = Sk_fg;
ck = Sk_fg(:,idx_col);
Ak_fg(:,idx_col) = [];

%% Build the matrix P
Pk = BuildP_BothDegree_SNTLN(m,m1,m2,n,n1,n2,k,k1,k2,alpha(ite),th1(ite),th2(ite),idx_col);
% Get the coefficients of f(x,y) in matrix form
% % fxy_vec = GetAsVector(fxy_matrix);
% % gxy_vec = GetAsVector(gxy_matrix);

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
x_ls = SolveAx_b(Ak_fg,ck);

% store this version of x_ls
% Insert a zero into the least squares solution x_ls so that 
% S_{k,k1,k2} x = c_{k,k1,k2}. Also so that when x is split into two 
% vectors x1 and x2. S(x1,x2) [f;g] =  ck.
first_part = x_ls(1:(idx_col-1));
second_part = x_ls(idx_col:end);
x = [first_part ; 0 ; second_part];

% Build Matrix Y, where Y(v,u)*[f;g] = S(f,g)*[u;v]
Yk = BuildY_BothDegree_SNTLN(x,m,m1,m2,n,n1,n2,k,k1,k2,alpha(ite),th1(ite),th2(ite));
% Test Y
%test1 = Y * [fxy_vec;gxy_vec];
%test2 = T * x;
%test1-test2;

% Calculate the initial residual r = ck - (Ak*x)
res_vec = ck - (Sk_fg*M*x_ls);

% Get the matrix p, which will store all the perturbations returned from LSE file
nEntries = nNonZeros_fxy...
    + nNonZeros_gxy ...
    + nCoeff_x ...
    + 3;


% Set the initial value of vector p to be zero
f = zeros(nNonZeros_fxy...
    + nNonZeros_gxy ...
    + nCoeff_x ...
    + 3,1);

%
% Set the intial value of E to the identity matrix
E = eye(nEntries);

%
% Create the matrix (T+N), initially N is empty so this is the same as T.
Tf = BuildT1_Both(fww_matrix,m,n-k,n1-k1,n2-k2);
Tg = BuildT1_Both(gww_matrix,n,m-k,m1-k1,m2-k2);
TN = [Tf alpha(ite).*Tg];

%
% Create The matrix (T+N) with respect to alpha
Tf_wrt_alpha = BuildT1_Both(Partial_fw_wrt_alpha_mat,m,n-k,n1-k1,n2-k2);
Tg_wrt_alpha = BuildT1_Both(Partial_alpha_gw_wrt_alpha_mat,n,m-k,m1-k1,m2-k2);
TN_wrt_alpha = [Tf_wrt_alpha Tg_wrt_alpha];

%
% Create The matrix (T+N) with respect to theta1
Tf_wrt_theta1 = BuildT1_Both(Partial_fw_wrt_theta1,m,n-k,n1-k1,n2-k2);
Tg_wrt_theta1 = BuildT1_Both(Partial_gw_wrt_theta1,n,m-k,m1-k1,m2-k2);
TN_wrt_theta1 = [Tf_wrt_theta1 alpha(ite) * Tg_wrt_theta1];

%
% Create The matrix (T+N) with respect to theta2
Tf_wrt_theta2 = BuildT1_Both(Partial_fw_wrt_theta2,m,n-k,n1-k1,n2-k2);
Tg_wrt_theta2 = BuildT1_Both(Partial_gw_wrt_theta2,n,m-k,m1-k1,m2-k2);
TN_wrt_theta2 = [Tf_wrt_theta2 alpha(ite) * Tg_wrt_theta2];

%%
% Create the matrix C for input into iteration

H_z     = Yk-Pk;

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
    
    % get the coefficients corresponding to f and g
    delta_zk        = y(1:nNonZeros_fxy + nNonZeros_gxy ,1);
    
    % Remove the zk coefficients from the list of coefficients
    y(1:nNonZeros_fxy + nNonZeros_gxy) = [];
    
    % Get the coefficients corresponding to x
    delta_xk        = y(1:nCoeff_x,1);
    
    % Remove them from the list of coefficients
    y(1:nCoeff_x) = [];
    
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
    Tf = BuildT1_Both(fww_matrix,m,n-k,n1-k1,n2-k2);
    Tg = BuildT1_Both(gww_matrix,n,m-k,m1-k1,m2-k2);
    Sk_fg = [Tf alpha(ite).*Tg];
    
    % Calculate the partial derivatives of fw and gw with respect to alpha
    Partial_fw_wrt_alpha_mat            = zeros(m1+1,m2+1);
    Partial_alpha_gw_wrt_alpha_mat      = gww_matrix;
    
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
    Tf_wrt_alpha = BuildT1_Both(Partial_fw_wrt_alpha_mat,m,n-k,n1-k1,n2-k2);
    Tg_wrt_alpha = BuildT1_Both(Partial_alpha_gw_wrt_alpha_mat,n,m-k,m1-k1,m2-k2);
    Partial_T_wrt_alpha = [Tf_wrt_alpha Tg_wrt_alpha];
    
    
    % Calculate the partial derivative of T with respect to theta1
    Tf_wrt_theta1 = BuildT1_Both(Partial_fw_wrt_theta1,m,n-k,n1-k1,n2-k2);
    Tg_wrt_theta1 = BuildT1_Both(Partial_gw_wrt_theta1,n,m-k,m1-k1,m2-k2);
    Partial_T_wrt_theta1 = [Tf_wrt_theta1 alpha(ite)*Tg_wrt_theta1];
    
    % Calculate the partial derivative of T with respect to theta2
    Tf_wrt_theta2 = BuildT1_Both(Partial_fw_wrt_theta2,m,n-k,n1-k1,n2-k2);
    Tg_wrt_theta2 = BuildT1_Both(Partial_gw_wrt_theta2,n,m-k,m1-k1,m2-k2);
    
    Partial_T_wrt_theta2 = [Tf_wrt_theta2 alpha(ite)*Tg_wrt_theta2];
    
    % Calculate the column c_{k} of DTQ that is moved to the right hand side
    ck = Sk_fg*e;
    
    % Calculate the derivatives of c_{k} with respect to \alpha and \theta
    Partial_ck_wrt_alpha        = Partial_T_wrt_alpha*e;
    Partial_ck_wrt_theta1       = Partial_T_wrt_theta1*e;
    Partial_ck_wrt_theta2       = Partial_T_wrt_theta2*e;
    
    % Create the vector of structured perturbations zf and zg applied
    % to F and G.
    z_fx      = zk(1:nNonZeros_fxy);
    z_gx      = zk(nNonZeros_fxy + 1 :end);
    
    % Calculate the entries of z_fw and z_gw
    
    z_fx_mat = GetAsMatrix([z_fx ; zeros(nZeros_fxy,1)],m1,m2);
    z_gx_mat = GetAsMatrix([z_gx ; zeros(nZeros_gxy,1)],n1,n2);
    
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
    N1 = BuildT1_Both(z_fw_mat,m,n-k,n1-k1,n2-k2);
    N2 = BuildT1_Both(z_gw_mat,n,m-k,m1-k1,m2-k2);
    N = [N1 alpha(ite).*N2];
    
    % Build the coefficient matrix N with respect to alpha
    Tf_wrt_alpha = BuildT1_Both(Partial_zfw_wrt_alpha,m,n-k,n1-k1,n2-k2);
    Tg_wrt_alpha = BuildT1_Both(Partial_zgw_wrt_alpha,n,m-k,m1-k1,m2-k2);
    Partial_N_wrt_alpha = [Tf_wrt_alpha Tg_wrt_alpha];
    
    
    % Calculate the derivatives of DNQ with respect to theta
    Tf_wrt_theta1 = BuildT1_Both(Partial_zfw_wrt_theta1,m,n-k,n1-k1,n2-k2);
    Tg_wrt_theta1 = BuildT1_Both(Partial_zgw_wrt_theta1,n,m-k,m1-k1,m2-k2);
    Partial_N_wrt_theta1 = [Tf_wrt_theta1 alpha(ite).*Tg_wrt_theta1];
    
    % Calculate the derivatives of DNQ with respect to theta
    Tf_wrt_theta2 = BuildT1_Both(Partial_zfw_wrt_theta2,m,n-k,n1-k1,n2-k2);
    Tg_wrt_theta2 = BuildT1_Both(Partial_zgw_wrt_theta2,n,m-k,m1-k1,m2-k2);
    Partial_N_wrt_theta2 = [Tf_wrt_theta2 alpha(ite).*Tg_wrt_theta2];
    
    % Calculate the column of DNQ that is moved to the right hand side, which
    % has the same structure as c_{k} the column of S_{k} moved to the RHS
    hk = N*e;
    
    % Calculate the derivative of h with respect to alpha
    h_alpha = Partial_N_wrt_alpha*e;
    
    % Calculate the derivative of h with respect to theta1
    h_theta1 = Partial_N_wrt_theta1*e;
    
    % Calculate the derivative of h with respect to theta2
    h_theta2 = Partial_N_wrt_theta2*e;
    
    % Build the matrix (T+N)
    
    Tf = BuildT1_Both(fww_matrix + z_fw_mat,m,n-k,n1-k1,n2-k2);
    Tg = BuildT1_Both(gww_matrix + z_gw_mat,n,m-k,m1-k1,m2-k2);
    TN = [Tf alpha(ite).*Tg];
    
    % Calculate the paritial derivative of (T+N) with respect to
    % alpha
    TN1_wrt_alpha = BuildT1_Both(Partial_fw_wrt_alpha_mat + Partial_zfw_wrt_alpha, m,n-k,n1-k1,n2-k2);
    TN2_wrt_alpha = BuildT1_Both(Partial_alpha_gw_wrt_alpha_mat + Partial_zgw_wrt_alpha,n,m-k, m1-k1,m2-k2);
    TN_alpha = [TN1_wrt_alpha TN2_wrt_alpha];
    
    
    % Calculate the paritial derivative of (T+N) with respect to
    % theta1
    TN1_wrt_theta1 = BuildT1_Both(Partial_fw_wrt_theta1 + Partial_zfw_wrt_theta1,m,n-k,n1-k1,n2-k2);
    TN2_wrt_theta1 = BuildT1_Both(Partial_gw_wrt_theta1 + Partial_zgw_wrt_theta1,n,m-k,m1-k1,m2-k2);
    TN_theta1 = [TN1_wrt_theta1 alpha(ite).*TN2_wrt_theta1];
    
    % Calculate the paritial derivative of (T+N) with respect to
    % theta2
    TN1_wrt_theta2 = BuildT1_Both(Partial_fw_wrt_theta2 + Partial_zfw_wrt_theta2,m,n-k,n1-k1,n2-k2);
    TN2_wrt_theta2 = BuildT1_Both(Partial_gw_wrt_theta2 + Partial_zgw_wrt_theta2,n,m-k,m1-k1,m2-k2);
    TN_theta2 = [TN1_wrt_theta2 alpha(ite).*TN2_wrt_theta2];
    
    % Update xk
    % Insert a zero into the least squares solution x_ls so that 
    % S_{k,k1,k2} x = c_{k,k1,k2}. Also so that when x is split into two 
    % vectors x1 and x2. S(x1,x2) [f;g] =  ck.
    first_part = xk(1:(idx_col-1));
    second_part = xk(idx_col:end);
    x = [first_part ; 0 ; second_part];

    % Calculate the matrix DY where Y is the Matrix such that E_{k}x = Y_{k}z.
    Yk = BuildY_BothDegree_SNTLN(x,m,m1,m2,n,n1,n2,k,k1,k2,alpha(ite),th1(ite),th2(ite));
    
    % Calculate the matrix P where ck = P * [f,g]
    Pk = BuildP_BothDegree_SNTLN(m,m1,m2,n,n1,n2,k,k1,k2,alpha(ite),th1(ite),th2(ite),idx_col);
    
    % Test P
    %ck;
    %ck2 = P*[fxy_vec;gxy_vec];
    %ck - ck2;
    
    
    %%
    
    % Get residual as a vector
    rk = (ck+hk) - TN*M*xk ;
    
    % Create the matrix C. This is made up of five submatrices, HZ, Hx,
    % H_alpha and H_theta1 and H_theta2.
    
    Hz          = Yk-Pk;
    
    Hx          = TN*M;
    
    H_alpha     = TN_alpha*M*xk - (Partial_ck_wrt_alpha + h_alpha);
    
    H_theta1    = TN_theta1*M*xk - (Partial_ck_wrt_theta1 + h_theta1);
    
    H_theta2    = TN_theta2*M*xk - (Partial_ck_wrt_theta2 + h_theta2);
    
    C = [Hz,Hx,H_alpha,H_theta1, H_theta2];  % the matrix C
    
    % Calculate the new right hand vector
    ek = ck + hk;
    
    % update gnew - used in LSE Problem.
    res_vec = rk;
    
    % Calculate the normalised residual of the solution.
    condition(ite) = norm(rk) / norm(ek);
    
    % Update fnew - used in LSE Problem.
    f = -(yy-start_point);
    
end

PlotGraphs_SNTLN();


%%

% Once iterations are complete, assign fx output, gx output, solution X
% output, alpha output and theta output.

% get the vector zk
zPert_f_vec = zk(1:nNonZeros_fxy);
zPert_f_mat = GetAsMatrix([zPert_f_vec ; zeros(nZeros_fxy,1)],m1,m2);

zPert_g_vec = zk(nNonZeros_fxy+1:end);
zPert_g_mat = GetAsMatrix([zPert_g_vec ; zeros(nZeros_gxy,1)],n1,n2);

% Set outputs of low rank approximation

fxy_lr = fxy_matrix + zPert_f_mat;

gxy_lr = gxy_matrix + zPert_g_mat;

X_lr  = xk;

alpha_lr = alpha(ite);

th1_lr = th1(ite);

th2_lr = th2(ite);

    


% Print the number of iterations
LineBreakLarge();
fprintf([mfilename sprintf('SNTLN converged within %i iterations \n', ite)]);
LineBreakLarge();

end

















