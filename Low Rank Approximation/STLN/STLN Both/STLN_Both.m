function [fxy_lr,gxy_lr] = STLN_Both(fxy_matrix,gxy_matrix,m,n,k,k1,k2,idx_col)
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
% m : total degree of f(x,y)
%
% n : total degree of g(x,y)
%
% k : Degree of d(x,y) 
%
% k1 : degree of d(x,y) with respect to x
%
% k2 : degree of d(x,y) with respect to y
%
% idx_col : index of optimal column for removal from S_{t_{1},t_{2}}(f,g)
%
% % Outputs 
%
% fxy_lr : Coefficients of f(x,y) with added perturbations
%
% gxy_lr : Coefficients of g(x,y) with added perturbations

global SETTINGS



% Get degree of polynomials f(x,y)
[m1,m2] = GetDegree(fxy_matrix);

% Get degree of polynomial g(x,y)
[n1,n2] = GetDegree(gxy_matrix);

% Get the number of coefficients in the polynomial f(x,y)
nCoeff_fxy = (m1+1) * (m2+1);
nNonZeros_fxy = GetNumNonZeros(m1,m2,m);
nZeros_fxy = nCoeff_fxy - nNonZeros_fxy;

% Get the number of coefficients in the polynomial g(x,y)
nCoeff_gxy = (n1+1) * (n2+1);
nNonZeros_gxy = GetNumNonZeros(n1,n2,n);
nZeros_gxy = nCoeff_gxy - nNonZeros_gxy;

% Get number of zeros in v
nNonZeros_vxy = GetNumNonZeros(n1-k1,n2-k2,n-k);

% Get number of zeros in u(x,y)
nNonZeros_uxy = GetNumNonZeros(m1-k1,m2-k2,m-k);

% Build the matrix T_{n1-k1,n2-k2}(f)
Tf = BuildT1_Both(fxy_matrix,m,n-k,n1-k1,n2-k2);

% Build the matrix T_{m1-k1,m2-k2}(g)
Tg = BuildT1_Both(gxy_matrix,n,m-k,m1-k1,m2-k2);

% Remove the columns of T(f) and T(g) which correspond to the zeros in u(x,y)
% and v(x,y) which are removed from the solution vector x.

% % 
% %
% Remove the extra rows of T1 and T2 associated with zeros of f*v and g*u
% Build the matrix S_{t,t_{1},t_{2}} with reduced rows and columns
St_fg = [Tf Tg];

% Remove optimal column
Ak_fg = St_fg;
Ak_fg(:,idx_col) = [];
ck = St_fg(:,idx_col);

% Build the matrix Et
Ak_zfzg = zeros(size(Ak_fg));
hk = zeros(size(ck));

% Build the vector z, which is the vector of perturbations of f and g.
z = zeros(nNonZeros_fxy + nNonZeros_gxy,1);

% Get the vector of coefficients zf
v_zf = z(1:nNonZeros_fxy);

% Get the vector of coefficeints zg
v_zg = z(nNonZeros_fxy + 1:end);

% Get zf as a matrix
% EDIT 10/03/2016 - vZ_fxy has zeros removed, so include the zeros to form
% a matrix m1+1 * m2+1
mat_zf = GetAsMatrix(...
    [
    v_zf; 
    zeros(nZeros_fxy,1)
    ],m1,m2);

% Get zg as a matrix
mat_zg = GetAsMatrix(...
    [
    v_zg; 
    zeros(nZeros_gxy,1)
    ],n1,n2);

% Get the vector x
% A_{t} x = c_{t}
xk = SolveAx_b(Ak_fg,ck);

x = ...
    [
    xk(1:idx_col-1);
    0;
    xk(idx_col:end);
    ];


% Build the matrix Y_{t}
% Where Y(x) * z = E(z) * x
Yk = BuildY_BothDegree_STLN(x,m,m1,m2,n,n1,n2,k,k1,k2);

% Test Y_{k}
v_fxy = GetAsVector(fxy_matrix);
v_fxy = v_fxy(1:nNonZeros_fxy,:);

v_gxy = GetAsVector(gxy_matrix);
v_gxy = v_gxy(1:nNonZeros_gxy,:);

test1 = Yk * [v_fxy;v_gxy];
test2 = Ak_fg * xk;
norm(test1-test2)

% Build the matrix P_{t}
% Where P * [f;g] = c_{t}
Pt = BuildP_BothDegree_STLN(m,m1,m2,n,n1,n2,k,k1,k2,idx_col);
test1 = Pt * [v_fxy;v_gxy];
test2 = ck;
norm(test1-test2)


% Get initial residual (A_{t}+E_{t})x = (c_{t} + h_{t})
xk = SolveAx_b(Ak_fg+Ak_zfzg,ck+hk);

H_z = Yk - Pt;

H_x = Ak_fg + Ak_zfzg;

C = [H_z H_x];

E = blkdiag( eye(nNonZeros_fxy + nNonZeros_gxy) , eye(nNonZeros_uxy + nNonZeros_vxy - 1));


% Get initial residual (A_{t}+E_{t})x = (c_{t} + h_{t})
res_vec = (ck + hk) - (Ak_fg*xk);

start_point     =   ...
    [...
    z;...
    xk;
    ];

yy = start_point;

f = -(yy-start_point);

% Initialise the iteration counter
ite = 1;

% Set the termination criterion
condition(ite) = norm(res_vec)./ norm(ck);

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
    xk    = xk + delta_xk;
    
    % Split vector z into vectors z_f and z_g
    v_zf = z(1:nNonZeros_fxy);
    v_zg = z(nNonZeros_fxy + 1:end);
    
    % Get zf as a matrix
    mat_zf = GetAsMatrix(...
        [
            v_zf; 
            zeros(nZeros_fxy,1)
        ]...
        ,m1,m2);
    
    % Get zg as a matrix
    mat_zg = GetAsMatrix(...
        [
            v_zg; 
            zeros(nZeros_gxy,1)
        ]...
        ,n1,n2);
    
    % Build the matrix B_{t} = [E1(zf) E2(zg)]
    T1_zf = BuildT1_Both(mat_zf,m,n-k,n1-k1,n2-k2);
    T1_zg = BuildT1_Both(mat_zg,n,m-k,m1-k1,m2-k2);
    
    % Build the matrix B_{t} equivalent to S_{t}
    St_zfzg = [T1_zf T1_zg];
    
    % Get the matrix E_{t} with optimal column removed
    Ak_zfzg = St_zfzg;
    Ak_zfzg(:,idx_col) = [];
    
    % Get the column vector h_{t}, the optimal column removed from B_{t},
    % and equivalent to c_{t} removed from S_{t}
    hk = St_zfzg(:,idx_col);
    
    % Get the updated vector x
    xk = SolveAx_b(Ak_fg + Ak_zfzg,ck + hk);
    
    x = [...
        xk(1:idx_col-1);...
        0;...
        xk(idx_col:end)];
    
    % Build the matrix Y_{t} where Y_{t}(x)*z = E_{t}(z) * x
    Yk = BuildY_BothDegree_STLN(x,m,m1,m2,n,n1,n2,k,k1,k2);
    
    % Get the residual vector
    res_vec = (ck + hk) - ((Ak_fg + Ak_zfzg) * xk);
  
    % Update the matrix C
    H_z = Yk - Pt;
    H_x = Ak_fg + Ak_zfzg;
    C = [H_z H_x];
    
    % Update fnew - used in LSE Problem.
    f = -(yy-start_point);
    
    % Update the termination criterion
    condition(ite) = norm(res_vec)./norm(ck+hk) ;
    
    
end

fprintf([mfilename ' : ' sprintf('\nRequired number of iterations: %i\n',ite)]);
    
fxy_lr = fxy_matrix + mat_zf;
gxy_lr = gxy_matrix + mat_zg;


PlotGraphs_STLN()


end

