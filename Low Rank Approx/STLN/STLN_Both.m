function [fxy_lr,gxy_lr] = STLN_Both(fxy_matrix,gxy_matrix,m,n,t,t1,t2,opt_col)
% Given coefficients f(x,y) and g(x,y) find the low rank approximation of
% the Syvlester subresultant S_{t_{1},t_{2}}.
%
% STLN(fxy,gxy,m,n,t1,t2,opt_col)
%
% % Inputs.
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
% % Outputs 
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

% Get the number of coefficients in v(x,y)
nCoeff_v = (n1-t1+1) * (n2-t2+1);

% Get the number of coefficients in u(x,y)
nCoeff_u = (m1-t1+1) * (m2-t2+1);

% Get the number of coefficients in the vector fv
nCoeff_fv = (m1+n1-t1+1) * (m2+n2-t2+1);

% Get the number of coefficients in the vector gu
nCoeff_gu = (n1+m1-t1+1) * (n2+m2-t2+1);

% Get number of zeros in v
nNonZeros_vxy = GetNumNonZeros(n1-t1,n2-t2,n-t);
nZeros_vxy = nCoeff_v - nNonZeros_vxy;

% Get number of zeros in u(x,y)
nNonZeros_uxy = GetNumNonZeros(m1-t1,m2-t2,m-t);
nZeros_uxy = nCoeff_u - nNonZeros_uxy;

% Get number of zeros in f(x,y)
nNonZeros_fxy = GetNumNonZeros(m1,m2,m);
nZeros_fxy = nCoeff_f - nNonZeros_fxy;

% Get number of zeros in g(x,y)
nNonZeros_gxy = GetNumNonZeros(n1,n2,n);
nZeros_gxy = nCoeff_g - nNonZeros_gxy;

% Get number of zeros in the product fv
nNonZeros_fv = GetNumNonZeros(m1+n1-t1,m2+n2-t2,m+n-t);
nZeros_fv = nCoeff_fv - nNonZeros_fv;

% Get number of zeros in the product gu
nNonZeros_gu = GetNumNonZeros(n1+m1-t1,n2+m2-t2,n+m-t);
nZeros_gu = nCoeff_gu - nNonZeros_gu;


% Get the number of coefficients in the unknown vector x, where A_{t}x =
% c_{t}.


% Build the Sylvester Matrix S(f,g)
T1 = BuildT1(fxy_matrix,n1-t1,n2-t2);
T2 = BuildT1(gxy_matrix,m1-t1,m2-t2);

% %
% %
% Remove the columns of T1 and T2 which correspond to the zeros in u(x,y)
% and v(x,y) which are removed from the solution vector x.
T1 = T1(:,1:nNonZeros_vxy);
T2 = T2(:,1:nNonZeros_uxy);

% % 
% %
% Remove the extra rows of T1 and T2 associated with zeros of f*v and g*u
T1 = T1(1:nNonZeros_fv,:);
T2 = T2(1:nNonZeros_gu,:);

% Build the matrix S_{t_{1},t_{2}} with reduced columns
St = [T1 T2];

% Remove optimal column
At = St;
At(:,opt_col) = [];
ct = St(:,opt_col);

% Build the matrix Et
Et = zeros(size(At));
ht = zeros(size(ct));

% EDIT 10/03/2016
% Build the vector z, which is the vector of perturbations of f and g.
z = zeros(nNonZeros_fxy + nNonZeros_gxy,1);

% Get the vector of coefficients zf
vZ_fxy = z(1:nNonZeros_fxy);

% Get the vector of coefficeints zg
vZ_gxy = z(nNonZeros_fxy + 1:end);

% Get zf as a matrix
% EDIT 10/03/2016 - vZ_fxy has zeros removed, so include the zeros to form
% a matrix m1+1 * m2+1
matZ_fxy = GetAsMatrix(...
    [
    vZ_fxy; 
    zeros(nZeros_fxy,1)
    ],m1,m2);

% Get zg as a matrix
matZ_gxy = GetAsMatrix(...
    [
    vZ_gxy; 
    zeros(nZeros_gxy,1)
    ],n1,n2);

% Get the vector x
% A_{t} x = c_{t}
x_ls = SolveAx_b(At,ct);

display(x_ls);


x = ...
    [
    x_ls(1:opt_col-1);
    0;
    x_ls(opt_col:end);
    ];


% Build the matrix Y_{t}
% Where Y(x) * z = E(z) * x
Yt = BuildYt(x,m,m1,m2,n,n1,n2,t,t1,t2);

v_fxy = GetAsVector(fxy_matrix);
v_fxy = v_fxy(1:nNonZeros_fxy,:);

v_gxy = GetAsVector(gxy_matrix);
v_gxy = v_gxy(1:nNonZeros_gxy,:);

display(Yt)
test1 = Yt * [v_fxy;v_gxy];
test2 = At * x_ls;
norm(test1-test2)

% Build the matrix P_{t}
% Where P * [f;g] = c_{t}
Pt = BuildPt(m,m1,m2,n,n1,n2,t,t1,t2,opt_col);
test1 = Pt * [v_fxy;v_gxy];
test2 = ct;
norm(test1-test2)


% Get initial residual (A_{t}+E_{t})x = (c_{t} + h_{t})
x_ls = SolveAx_b(At+Et,ct+ht);

H_z = Yt - Pt;

H_x = At + Et;

C = [H_z H_x];

E = blkdiag( eye(nNonZeros_fxy + nNonZeros_gxy) , zeros(nNonZeros_uxy + nNonZeros_vxy - 1));


% Get initial residual (A_{t}+E_{t})x = (c_{t} + h_{t})
res_vec = (ct + ht) - (At*x_ls);

start_point     =   ...
    [...
    z;...
    x_ls;
    ];

yy = start_point;

f = -(yy-start_point);

% Initialise the iteration counter
ite = 1;

% Set the termination criterion
condition(ite) = norm(res_vec)./ norm(ct);

while condition(ite) >  SETTINGS.MAX_ERROR_SNTLN &&  ite < SETTINGS.MAX_ITERATIONS_SNTLN
    
    % Increment interation counter
    ite = ite + 1;
    
    % Get small petrubations by LSE problem ||Ey-f||  subject to  Cy=g
    y_lse = LSE(E,f,C,res_vec);
    
    % Increment cummulative peturbations
    yy = yy + y_lse;
    
    % obtain the small changes to z and x
    nEntries_z      = nNonZeros_fxy + nNonZeros_gxy;
    delta_zk        = y_lse(1:nEntries_z);
    delta_xk        = y_lse((nEntries_z+1):end);
    
    % Update z and x
    z       = z + delta_zk;
    x_ls    = x_ls + delta_xk;
    
    % Split vector z into vectors z_f and z_g
    vZ_fxy = z(1:nNonZeros_fxy);
    vZ_gxy = z(nNonZeros_fxy + 1:end);
    
    % Get zf as a matrix
    matZ_fxy = GetAsMatrix(...
        [
            vZ_fxy; 
            zeros(nZeros_fxy,1)
        ]...
        ,m1,m2);
    
    % Get zg as a matrix
    matZ_gxy = GetAsMatrix(...
        [
            vZ_gxy; 
            zeros(nZeros_gxy,1)
        ]...
        ,n1,n2);
    
    % Build the matrix B_{t} = [E1(zf) E2(zg)]
    E1 = BuildT1(matZ_fxy,n1-t1,n2-t2);
    E2 = BuildT1(matZ_gxy,m1-t1,m2-t2);
    
    % Remove the extra columns of E1 and E2
    E1 = E1(:,1:nNonZeros_vxy);
    E2 = E2(:,1:nNonZeros_uxy);
    
    % Remove Extra rows of E1 and E2
    E1 = E1(1:nNonZeros_fv,:);
    E2 = E2(1:nNonZeros_gu,:);
    
    % Build the matrix B_{t} equivalent to S_{t}
    Bt = [E1 E2];
    
    % Get the matrix E_{t} with optimal column removed
    Et = Bt;
    Et(:,opt_col) = [];
    
    % Get the column vector h_{t}, the optimal column removed from B_{t},
    % and equivalent to c_{t} removed from S_{t}
    ht = Bt(:,opt_col);
    
    % Get the updated vector x
    x_ls = SolveAx_b(At + Et,ct + ht);
    
    x = [...
        x_ls(1:opt_col-1);...
        0;...
        x_ls(opt_col:end)];
    
    % Build the matrix Y_{t} where Y_{t}(x)*z = E_{t}(z) * x
    Yt = BuildYt(x,m,m1,m2,n,n1,n2,t,t1,t2);
    
    % Get the residual vector
    res_vec = (ct+ht) - ((At+Et)*x_ls);
  
    % Update the matrix C
    H_z = Yt - Pt;
    H_x = At + Et;
    C = [H_z H_x];
    
    % Update fnew - used in LSE Problem.
    f = -(yy-start_point);
    
    % Update the termination criterion
    condition(ite) = norm(res_vec)./norm(ct+ht) ;
    
    
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
    
    fxy_lr = fxy_matrix + matZ_fxy;
    gxy_lr = gxy_matrix + matZ_gxy;
else
    fxy_lr = fxy_matrix;
    gxy_lr = gxy_matrix;
end


end

function Yt = BuildYt(x,m,m1,m2,n,n1,n2,t,t1,t2)
% Build the matrix Y_{t}
% Where Y(x) * z = E(z) * x
% The vector x only contains the non-zero entries of v(x,y) and u(x,y)

nCoefficients_uxy = (m1-t1+1) * (m2-t2+1);
nCoefficients_vxy = (n1-t1+1) * (n2-t2+1);

nNonZeros_uxy = GetNumNonZeros(m1-t1,m2-t2,m-t);
nNonZeros_vxy = GetNumNonZeros(n1-t1,n2-t2,n-t);

nZeros_uxy = nCoefficients_uxy - nNonZeros_uxy;
nZeros_vxy = nCoefficients_vxy - nNonZeros_vxy;

nNonZeros_fxy = GetNumNonZeros(m1,m2,m);
nNonZeros_gxy = GetNumNonZeros(n1,n2,n);

xv = x(1:nNonZeros_vxy);
xu = 1.* x(nNonZeros_vxy + 1 : end);

% get x_u as a matrix
mat_xu = GetAsMatrix(...
    [
        xu;
        zeros(nZeros_uxy,1)
    ]...
    ,m1-t1,m2-t2);

% Get x_v as a matrix
mat_xv = GetAsMatrix(...
    [
        xv;
        zeros(nZeros_vxy,1);
    ]...
    ,n1-t1,n2-t2);

% Build the matrix C1 and C2
C1 = BuildT1(mat_xv,m1,m2);
C2 = BuildT1(mat_xu,n1,n2);

% Remove the columns of C1(v) and C2(u) corresponding to zeros in f(x,y) and
% g(x,y)
C1 = C1(:,1:nNonZeros_fxy);
C2 = C2(:,1:nNonZeros_gxy);

% Remove the rows of C1(v) and C2(u) corresponding to zeros in the product
% f(x,y)*v(x,y) and g(x,y)*u(x,y)

nNonZeros_fv = GetNumNonZeros(m1+n1-t1,m2+n2-t2,m+n-t);
nNonZeros_gu = GetNumNonZeros(n1+m1-t1,n2+m2-t2,n+m-t);

C1 = C1(1:nNonZeros_fv,:);
C2 = C2(1:nNonZeros_gu,:);

Yt = [C1 C2];



end

function Pt = BuildPt(m,m1,m2,n,n1,n2,t,t1,t2,opt_col_index)
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
% t
%
% t1 :
%
% t2 :


% Get the number of coefficients in polynomial f(x,y)
nCoeff_f = (m1+1).*(m2+1);

% Get the number of coefficients in polynomial g(x,y)
nCoeff_g = (n1+1).*(n2+1);

% Get the number of nonZero coefficients in f(x,y)
nNonZeros_fxy = GetNumNonZeros(m1,m2,m);

% Get the number of nonZeros coefficients in g(x,y)
nNonZeros_gxy = GetNumNonZeros(n1,n2,n);

% Number of columns in T1 of the sylvester matrix
nColumnsT1 = nNonZeros_fxy;



if opt_col_index <= nColumnsT1
    % Optimal column in first partition
    
    % % Build the matrix P
    
    % Build the matrix P1
    P1 = BuildPt_sub(m,m1,m2,n,n1,n2,t,t1,t2,opt_col_index);
    
    % Build the matrix P2
    nRows = (m1+n1-t1+1)*(m2+n2-t2+1);
    nNonZeros_fv = GetNumNonZeros(m1+n1-t1,m2+n2-t2,m+n-t);
    P2 = zeros(nNonZeros_fv,nNonZeros_gxy);
    
    
    
    
else
    % Optimal column in second partition
    
    % Build the matrix P1
    nRows = (m1+n1-t1+1)*(m2+n2-t2+1);
    nNonZeros_fv = GetNumNonZeros(m1+n1-t1,m2+n2-t2,m+n-t);
    
    P1 = zeros(nNonZeros_fv,nNonZeros_fxy);
    
    
    
    % Get the position of the optimal column with respect to T(g), so
    % remove the number of columns in T(f).
    opt_col_rel = opt_col_index - nColumnsT1;
    
    % Build the matrix P2
    P2 = BuildPt_sub(n,n1,n2,m,m1,m2,t,t1,t2,opt_col_rel);
    
    
    
    % Build the matrix P.
    
    
end

Pt = [P1 P2];

end

function P = BuildPt_sub(m,m1,m2,n,n1,n2,t,t1,t2,opt_col)
% BuildPt_sub(m1,m2,n1,n2,opt_col,t1,t2)
%
% Build the matrix P, used in STLN function. P is a matrix which is
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


% Get number of nonzeros of f(x,y)
nNonZeros_fxy = GetNumNonZeros(m1,m2,m);

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

%
ihat = i+1;
jhat = j+1;

%
nRows_f = m1+1;
nCols_f = m2+1;

% insert the theta matrix into the zero matrix
padd_mat(ihat:i+nRows_f, jhat:j+nCols_f) = mat;

% Get the padded matrix as a vector
vec_padd_mat = GetAsVector(padd_mat);

% Diagonalise the vector.
diag_mat_vec_padd_mat = diag(vec_padd_mat);

% Remove the zero columns
diag_mat_vec_padd_mat( :, ~any(diag_mat_vec_padd_mat,1) ) = [];  %columns


P = diag_mat_vec_padd_mat;

% Remove the columns corresponding to zeros in f
P = P(:,1:nNonZeros_fxy);

% Remove the rows corresponding to zeros of f*v
nNonZeros_fv = GetNumNonZeros(m1+n1-t1,m2+n2-t2,m+n-t);

P = P(1:nNonZeros_fv,:);


end