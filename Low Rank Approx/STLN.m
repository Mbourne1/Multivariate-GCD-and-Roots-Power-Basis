function [fxy_matrix_out,gxy_matrix_out] = STLN(fxy_matrix,gxy_matrix,m,n,t1,t2,opt_col)
% Given coefficients f(x,y) and g(x,y) find the low rank approximation of
% the Syvlester subresultant S_{t_{1},t_{2}}.
%
% STLN(fxy,gxy,m,n,t1,t2,opt_col)
%
% Inputs.
%
% fxy_matrix : Coefficients of polynomial f(x,y)
%
% gxy_matrix : Coefficients of polynomial g(x,y)
%
% m : total degree of f
%
% n : total degree of g
%
% t1 : degree of d(x,y) with respect to x
%
% t2 : degree of d(x,y) with respect to y
%
% opt_col : index of optimal column for removal from S_{t_{1},t_{2}}(f,g)
%
% Outputs 
%
% fxy_matrix : Coefficients of f(x,y) with added perturbations
global SETTINGS


% Get degree of polynomials f(x,y)
[m1,m2] = GetDegree(fxy_matrix);

% Get degree of polynomial g(x,y)
[n1,n2] = GetDegree(gxy_matrix);

% Get the number of coefficients in the polynomial f(x,y)
nCoeff_f = (m1+1) * (m2+1);

% Get the number of coefficients in the polynomial g(x,y)
nCoeff_g = (n1+1) * (n2+1);

% Get the number of coefficients in both f(x,y) and g(x,y)
%num_coeff = num_coeff_f + num_coeff_g;

% Get the number of coefficients in v(x,y)
nCoeff_v = (n1-t1+1) * (n2-t2+1);

% Get the number of coefficients in u(x,y)
nCoeff_u = (m1-t1+1) * (m2-t2+1);

% Since we know m,m1,m2,n,n1,n2 we can remove columns corresponding to the
% coefficients which will definitely be zeros.
SETTINGS.BOOL_REMOVECOLS = 'n';

switch SETTINGS.BOOL_REMOVECOLS
    case 'y'
        % % Get number of zeros in f(x,y)
        diff_f = m1 + m2 - m;
        diff_g = n1 + n2 - n;
        
        if diff_f ~=0
            nZeros_f = nchoosek(diff_f + 1,2);
        else
            nZeros_f = 0;
        end
        
        if diff_g ~=0
            nZeros_g = nchoosek(diff_g + 1,2);
        else
            nZeros_g = 0;
        end
    case 'n'
        
        nZeros_f = 0 ;
        nZeros_g = 0;
end

% Get the number of coefficients in the unknown vector x, where A_{t}x =
% c_{t}.
%num_coeff_x = num_coeff_u + num_coeff_v - 1;

% Build the Sylvester Matrix S(f,g)
T1 = BuildT1(fxy_matrix,n1-t1,n2-t2);
T2 = BuildT1(gxy_matrix,m1-t1,m2-t2);
St = [T1 T2];

% Remove optimal column
At = St;
At(:,opt_col) = [];
ct = St(:,opt_col);

% Build the matrix Et
Et = zeros(size(At));
ht = zeros(size(ct));

% EDIT 10/03/2016
z = zeros(nCoeff_f + nCoeff_g - nZeros_f - nZeros_g,1);

% Get the vector of coefficients zf
vZ_fxy = z(1:nCoeff_f - nZeros_f);

% Get the vector of coefficeints zg
vZ_gxy = z(nCoeff_f-nZeros_f+1:end);

% Get zf as a matrix
% EDIT 10/03/2016 - vZ_fxy has zeros removed, so include the zeros to form
% a matrix m1+1 * m2+1
matZ_fxy = GetAsMatrix([vZ_fxy; zeros(nZeros_f,1)],m1,m2);

% Get zg as a matrix
matZ_gxy = GetAsMatrix([vZ_gxy; zeros(nZeros_g,1)],n1,n2);

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
Yt = BuildYt(x,m,m1,m2,n,n1,n2,t1,t2);
v_fxy = GetAsVector(fxy_matrix);
v_fxy(end-(nZeros_f-1):end) = [];

v_gxy = GetAsVector(gxy_matrix);
v_gxy(end-(nZeros_g-1):end) = [];

test1 = Yt * [v_fxy;v_gxy];
test2 = At * x_ls;

norm(test1-test2)

% Build the matrix P_{t}
% Where P * [f;g] = c_{t}
Pt = BuildPt(m,m1,m2,n,n1,n2,opt_col,t1,t2);

test1 = Pt * [v_fxy;v_gxy];
test2 = ct;

norm(test1-test2)


% Get initial residual (A_{t}+E_{t})x = (c_{t} + h_{t})
x_ls = SolveAx_b(At+Et,ct+ht);

% Get the residual vector
g = (ct + ht) - (At+Et)*x_ls;

H_z = Yt - Pt;

H_x = At + Et;

C = [H_z H_x];

% Build the matrix E
nEntries = (nCoeff_f-nZeros_f) + (nCoeff_g-nZeros_g) + nCoeff_u + nCoeff_v - 1;

E = blkdiag( eye(nCoeff_f -  + nCoeff_g) , zeros(nCoeff_u + nCoeff_v - 1));
%E = eye(nEntries);

yy = zeros((nCoeff_f-nZeros_f) + (nCoeff_g-nZeros_g) + nCoeff_u + nCoeff_v - 1,1);

% Get initial residual (A_{t}+E_{t})x = (c_{t} + h_{t})
g = (ct + ht) - At*x_ls;

start_point     =   ...
    [...
    z;...
    x_ls;
    ];

f = -(yy-start_point);


% Initialise the iteration counter
ite = 1;

% Set the termination criterion
condition(ite) = norm(g)./ norm(ct);

while condition(ite) >  SETTINGS.MAX_ERROR_SNTLN &&  ite < SETTINGS.MAX_ITERATIONS_SNTLN
    
    % Increment interation counter
    ite = ite + 1;
    
    % Get small petrubations by LSE problem ||Ey-f||  subject to  Cy=g
    y_lse = LSE(E,f,C,g);
    
    % Increment cummulative peturbations
    yy = yy + y_lse;
    
    % obtain the small changes to z and x
    nEntries_z      = nCoeff_f + nCoeff_g - nZeros_f - nZeros_g;
    delta_zk        = y_lse(1:nEntries_z);
    delta_xk        = y_lse((nEntries_z+1):end);
    
    % Update z and x
    z       = z + delta_zk;
    x_ls    = x_ls + delta_xk;
    
    % Split vector z into vectors z_f and z_g
    vZ_fxy = z(1:nCoeff_f-nZeros_f);
    vZ_gxy = z(nCoeff_f-nZeros_f + 1:end);
    
    % Get zf as a matrix
    matZ_fxy = GetAsMatrix([vZ_fxy; zeros(nZeros_f,1)],m1,m2);
    
    % Get zg as a matrix
    matZ_gxy = GetAsMatrix([vZ_gxy; zeros(nZeros_g,1)],n1,n2);
    
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
    %x_ls = SolveAx_b(At+Et,ct+ht);
    
    x = [...
        x_ls(1:opt_col-1);...
        0;...
        x_ls(opt_col:end)];
    
    % Build the matrix Y_{t} where Y_{t}(x)*z = E_{t}(z) * x
    Yt = BuildYt(x,m,m1,m2,n,n1,n2,t1,t2);
    
    % Get the residual vector
    g = (ct+ht) - ((At+Et)*x_ls);
    
    % Update the matrix C
    H_z = Yt - Pt;
    H_x = At + Et;
    C = [H_z H_x];
    
    % Update fnew - used in LSE Problem.
    f = -(yy-start_point);
    
    % Update the termination criterion
    condition(ite) = norm(g)./norm(ct+ht) ;
    
    
end

fprintf('\nRequired number of iterations: %i\n',ite)

switch SETTINGS.PLOT_GRAPHS
    case 'y'
        
        title = sprintf('%s - Residuals',mfilename());
        figure('name',title)
        hold on
        plot(log10(condition),'-s')
        hold off
    case 'n'
end


if condition(ite) < condition(1)
    
    fxy_matrix_out = fxy_matrix + matZ_fxy;
    gxy_matrix_out = gxy_matrix + matZ_gxy;
else
    fxy_matrix_out = fxy_matrix;
    gxy_matrix_out = gxy_matrix;
end

display([GetAsVector(fxy_matrix) GetAsVector(matZ_fxy)]);
display([GetAsVector(gxy_matrix) GetAsVector(matZ_gxy)]);

end

function Yt = BuildYt(x,m,m1,m2,n,n1,n2,t1,t2)
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

% Remove the columns of C1 which correspond to zero values in the vector
% of coefficients of f

diff_f = m1 + m2 - m;
diff_g = n1 + n2 - n;

global SETTINGS
switch SETTINGS.BOOL_REMOVECOLS
    case 'y'
        
        if diff_f ~=0
            num_zeros_f = nchoosek(diff_f + 1,2);
        else
            num_zeros_f = 0;
        end
        
        if diff_g ~=0
            num_zeros_g = nchoosek(diff_g + 1,2);
        else
            num_zeros_g = 0;
        end
        
    case 'n'
        num_zeros_f = 0;
        num_zeros_g = 0;
    otherwise
        error('err')
end
C1(:,end-(num_zeros_f-1) : end) = [];
C2(:,end-(num_zeros_g-1) : end) = [];

Yt = [C1 C2];



end

function Pt = BuildPt(m,m1,m2,n,n1,n2,opt_col_index,t1,t2)
% BuildPt(m,m1,m2,n,n1,n2,opt_col,t1,t2)
%
% Build the matrix P_{t}, such that the matrix vector product P*[f;g] gives
% the column c_{t}.
%
% P_{t} * [f;g] = c_{t}
%
% Inputs
%
% m  :
%
% m1 :
%
% m2 :
%
% n  :
%
% n1 :
%
% n2 :
%
% opt_col :
%
% t1 :
%
% t2 :


% Get the number of coefficients in polynomial f
nCoeff_f = (m1+1).*(m2+1);

% Get the number of coefficients in polynomial g
nCoeff_g = (n1+1).*(n2+1);

% Number of columns in T1 of the sylvester matrix
nColumnsT1 = (n1-t1+1) * (n2-t2+1);

if opt_col_index <= nColumnsT1
    % Optimal column in first partition
    
    % % Build the matrix P
    
    % Build the matrix P1
    P1 = BuildPt_sub(m1,m2,n1,n2,opt_col_index,t1,t2);
    
    % Build the matrix P2
    rows = (m1+n1-t1+1)*(m2+n2-t2+1);
    P2 = zeros(rows,nCoeff_g);
    
    
    
else
    % Optimal column in second partition
    
    % Build the matrix P1
    rows = (m1+n1-t1+1)*(m2+n2-t2+1);
    P1 = zeros(rows,nCoeff_f);
    
    % Build the matrix P2
    % Get the position of the optimal column with respect to T(g)
    opt_col_rel = opt_col_index - nColumnsT1;
    P2 = BuildPt_sub(n1,n2,m1,m2,opt_col_rel,t1,t2);
    
    % Build the matrix P.
    
    
end
diff_f = m1 + m2 - m;
diff_g = n1 + n2 - n;

global SETTINGS
switch SETTINGS.BOOL_REMOVECOLS
    case 'y'
        if diff_f ~=0
            num_zeros_f = nchoosek(diff_f + 1,2);
        else
            num_zeros_f = 0;
        end
        
        if diff_g ~=0
            num_zeros_g = nchoosek(diff_g + 1,2);
        else
            num_zeros_g = 0;
        end
    case 'n'
        num_zeros_f = 0;
        num_zeros_g = 0;
end

P1(:,end-(num_zeros_f-1) : end) = [];
P2(:,end-(num_zeros_g-1) : end) = [];

Pt = [P1 P2];

end

function P = BuildPt_sub(m1,m2,n1,n2,opt_col,t1,t2)
% BuildPt_sub(m1,m2,n1,n2,opt_col,t1,t2)
%
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

nRows_f = m1+1;
nCols_f = m2+1;

% inser the theta matrix into the zero matrix
padd_mat(ihat:i+nRows_f, jhat:j+nCols_f) = mat;

% Get the padded matrix as a vector
vec_padd_mat = GetAsVector(padd_mat);

% Diagonalise the vector.
diag_mat_vec_padd_mat = diag(vec_padd_mat);


% Remove the zero columns
diag_mat_vec_padd_mat( :, ~any(diag_mat_vec_padd_mat,1) ) = [];  %columns


P = diag_mat_vec_padd_mat;

end