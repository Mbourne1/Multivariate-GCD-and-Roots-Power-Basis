function [fxy,gxy] = STLN(fxy,gxy,t1,t2,opt_col)
% Given coefficients f(x,y) and g(x,y) find the low rank approximation of
% the Syvlester subresultant S_{t_{1},t_{2}}.

global MAX_ERROR_SNTLN
global MAX_ITERATIONS_SNTLN
global PLOT_GRAPHS

% Get degree of polynomials f(x,y)
[r,c] = size(fxy);
m1 = r - 1;
m2 = c - 1;

% Get degree of polynomial g(x,y)
[r,c] = size(gxy);
n1 = r - 1;
n2 = c - 1;

% Get the number of coefficients in the polynomial f(x,y)
num_coeff_f = (m1+1) * (m2+1);

% Get the number of coefficients in the polynomial g(x,y)
num_coeff_g = (n1+1) * (n2+1);

% Get the number of coefficients in both f(x,y) and g(x,y)
%num_coeff = num_coeff_f + num_coeff_g;

% Get the number of coefficients in v(x,y)
num_coeff_v = (n1-t1+1) * (n2-t2+1);

% Get the number of coefficients in u(x,y)
num_coeff_u = (m1-t1+1) * (m2-t2+1);

% Get the number of coefficients in the unknown vector x, where A_{t}x =
% c_{t}.
%num_coeff_x = num_coeff_u + num_coeff_v - 1;

% Build the Sylvester Matrix S(f,g)
T1 = BuildT1(fxy,n1-t1,n2-t2);
T2 = BuildT1(gxy,m1-t1,m2-t2);
St = [T1 T2];

% Remove optimal column
At = St;
At(:,opt_col) = [];
ct = St(:,opt_col);

% Build the matrix Et
Et = zeros(size(At));
ht = zeros(size(ct));

z = zeros(num_coeff_f + num_coeff_g,1);
% Get the vector of coefficients zf
vZ_fxy = z(1:num_coeff_f);
% Get the vector of coefficeints zg
vZ_gxy = z(num_coeff_f+1:end);

% Get zf as a matrix
matZ_fxy = GetAsMatrix(vZ_fxy,m1,m2);
% Get zg as a matrix
matZ_gxy = GetAsMatrix(vZ_gxy,n1,n2);

% Get the vector x
% A_{t} x = c_{t}
x_ls = SolveAx_b(At,ct);

x = ...
    [
    x_ls(1:opt_col-1);
    0;
    x_ls(opt_col:end);
    ];


% Build the matrix Y_{t}
% Where Y(x) * z = E(z) * x

Yt = BuildY_RelativeDegree(x,m1,m2,n1,n2,t1,t2);
%v_fxy = GetAsVector(fxy);
%v_gxy = GetAsVector(gxy);

%test1 = Yt * [v_fxy;v_gxy];
%test2 = At * x_ls;

% Build the matrix P_{t}
% Where P * [f;g] = c_{t}
Pt = BuildP_RelativeDegree(m1,m2,n1,n2,opt_col,t1,t2);

%test1 = Pt * [v_fxy;v_gxy];
%test2 = ct;

% Get initial residual (A_{t}+E_{t})x = (c_{t} + h_{t})
x_ls = pinv(At+Et) * (ct + ht);
g = (ct + ht) - (At+Et)*x_ls;


H_z = Yt - Pt;

H_x = At + Et;

C = [H_z H_x];

% Build the matrix E
%nEntries = num_coeff_f + num_coeff_g + num_coeff_u + num_coeff_v - 1;


E = blkdiag( eye(num_coeff_f + num_coeff_g) , zeros(num_coeff_u + num_coeff_v - 1));
%E = eye(nEntries);

yy = zeros(num_coeff_f + num_coeff_g + num_coeff_u + num_coeff_v - 1,1);

start_point     =   ...
    [...
    z;...
    x_ls;
    ];

f = -(yy - start_point);


% Initialise the iteration counter
ite = 1;

% Set the termination criterion
condition(ite) = norm(g);

while condition(ite) >  MAX_ERROR_SNTLN &&  ite < MAX_ITERATIONS_SNTLN
    
    % Increment interation counter
    ite = ite + 1;
    
    % Get small petrubations by LSE
    y_lse = LSE(E,f,C,g);
    
    % Increment cummulative peturbations
    yy = yy + y_lse;
    
    % obtain the small changes to z and x
    nEntries_z      = num_coeff_f + num_coeff_g;
    delta_zk        = y_lse(1:nEntries_z);
    delta_xk        = y_lse((nEntries_z+1):end);
    
    % Update z and x
    z       = z + delta_zk;
    %x_ls    = x_ls + delta_xk;
    
    % Split vector z into vectors z_f and z_g
    vZ_fxy = z(1:num_coeff_f);
    vZ_gxy = z(num_coeff_f + 1:end);
    
    % Get zf as a matrix
    matZ_fxy = GetAsMatrix(vZ_fxy,m1,m2);
    
    % Get zg as a matrix
    matZ_gxy = GetAsMatrix(vZ_gxy,n1,n2);
    
    % Build the matrix B_{t} = [E1(zf) E2(zg)]
    E1 = BuildT1(matZ_fxy,n1-t1,n2-t2);
    E2 = BuildT1(matZ_gxy,m1-t1,m2-t2);
    
    % Build the matrix B_{t} equivalent to S_{t}
    Bt = [E1 E2];
    
    % Get the matrix E_{t} with optimal column removed
    Et = Bt;
    Et(:,opt_col) = [];
    
    % Get the column vector h_{t}, the optimal column removed from B_{t},
    % and equivalent to c_{t} removed from S_{t}
    ht = Bt(:,opt_col);
    
    % Get the updated vector x
    x_ls = SolveAx_b(At+Et,ct+ht);
    x = [x_ls(1:opt_col-1) ; 0 ; x_ls(opt_col:end)];
    
    % Build the matrix Y_{t} where Y_{t}(x)*z = E_{t}(z) * x
    Yt = BuildYt(x,m1,m2,n1,n2,t1,t2);
    
    % Get the residual vector
    g = (ct+ht) - ((At+Et)*x_ls);
    
    % Update the matrix C
    H_z = Yt - Pt;
    H_x = At + Et;
    C = [H_z H_x];
    
    % Update f - used in LSE Problem.
    f = -(yy+start_point);
    
    % Update the termination criterion
    condition(ite) = norm(g) ;
    
    
end

fprintf('\nRequired number of iterations: %i\n',ite)

switch PLOT_GRAPHS
    case 'y'
        
        figure('name','STLN - Residuals')
        hold on
        plot(log10(condition),'-s')
        hold off
    case 'n'
end


if cond(ite) < cond(1)
    
    fxy = fxy + matZ_fxy;
    gxy = gxy + matZ_gxy;
else
    %fxy = fxy;
    %gxy = gxy;
end


end

function Yt = BuildY_RelativeDegree(x,m1,m2,n1,n2,t1,t2)
% Build the matrix Y_{t}
% Where Y(x) * z = E(z) * x

%nCoefficients_xu = (m1-t1+1) * (m2-t2+1);
nCoefficients_xv = (n1-t1+1) * (n2-t2+1);

xv = x(1:nCoefficients_xv);
xu = x(nCoefficients_xv+1:end);

% get x_u as a matrix
mat_xu = GetAsMatrix(xu,m1-t1,m2-t2);
mat_xv = GetAsMatrix(xv,n1-t1,n2-t2);

C1 = BuildT1(mat_xv,m1,m2);
C2 = BuildT1(mat_xu,n1,n2);

Yt = [C1 C2];

end

function Pt = BuildP_RelativeDegree(m1,m2,n1,n2,opt_col,t1,t2)
% Build the matrix P_{t}
% Where P * [f;g] = c_{t}


% Get the number of coefficients in polynomial f
num_coeff_f = (m1+1).*(m2+1);

% Get the number of coefficients in polynomial g
num_coeff_g = (n1+1).*(n2+1);

% Number of columns in T1 of the sylvester matrix
num_cols_T1 = (n1-t1+1) * (n2-t2+1);

if opt_col <= num_cols_T1
    % Optimal column in first partition
    
    % % Build the matrix P
    
    % Build the matrix P1
    P1 = BuildP1_RelativeDegree(m1,m2,n1,n2,opt_col,t1,t2);
    
    % Build the matrix P2
    rows = (m1+n1-t1+1)*(m2+n2-t2+1);
    P2 = zeros(rows,num_coeff_g);
    
    % Build the matrix P
    Pt = [P1 P2];
    
else
    % Optimal column in second partition
    
    % Build the matrix P1
    rows = (m1+n1-t1+1)*(m2+n2-t2+1);
    P1 = zeros(rows,num_coeff_f);
    
    % Build the matrix P2
    % Get the position of the optimal column with respect to T(g)
    opt_col_rel = opt_col - num_cols_T1;
    P2 = BuildP1_RelativeDegree(n1,n2,m1,m2,opt_col_rel,t1,t2);
    
    % Build the matrix P.
    Pt = [P1 P2];
    
end

end

function P = BuildP1_RelativeDegree(m1,m2,n1,n2,opt_col,t1,t2)
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
mat = ones(m1+1,m2+1);


% pad with zeros
%num_mult_wrt_x = n1-t1;
%num_mult_wrt_y = n2-t2;

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