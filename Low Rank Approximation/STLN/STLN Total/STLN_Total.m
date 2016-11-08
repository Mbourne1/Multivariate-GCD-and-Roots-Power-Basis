function [fxy_lr,gxy_lr] = STLN_Total(fxy_matrix,gxy_matrix,m,n,k,idx_col)
% STLN(fxy_matrix, gxy_matrix, m, n, k, idx_col)
%
% Given coefficients f(x,y) and g(x,y) find the low rank approximation of
% the Syvlester subresultant S_{k}.
%
%
% Inputs.
%
% fxy_matrix : Coefficients of polynomial f(x,y)
%
% gxy_matrix : Coefficients of polynomial g(x,y)
%
% m : Total degree of f(x,y)
%
% n : Total degree of g(x,y)
%
% k : Total degree of d(x,y) 
%
% idx_col : Index of optimal column for removal from S_{k_{1},k_{2}}(f,g)
%
% Outputs
%
% fxy_lr : Coefficients of f(x,y) with added perturbations
%
% gxy_lr : Coefficients of g(x,y) with added perturbations

%%
% pad the coefficients of fxy and gxy
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
global SETTINGS

% Get the number of coefficients in the polynomial f(x,y)
nNonZeros_fxy = nchoosek(m+2,2);
nZeros_fxy = nchoosek(m+1,2);

% Get the number of coefficients in the polynomial g(x,y)
nNonZeros_gxy = nchoosek(n+2,2);
nZeros_gxy = nchoosek(n+1,2);

% Get the number of coefficients in v(x,y)
nNonZeros_vxy = nchoosek(n-k+2,2);

% Get the number of coefficients in u(x,y)
nNonZeros_uxy = nchoosek(m-k+2,2);

% Build the Sylvester Matrix S_{k1,k2}(f,g)

% Build the partiton T_{n1-k1,n2-k2}(f)
Tk_fxy = BuildT1_Total(fxy_matrix,m,n-k);

% Build the partition T_{m1-k1,m2-k2}(g)
Tk_gxy = BuildT1_Total(gxy_matrix,n,m-k);

% Build the Sylvester subresultant matrix S_{k1,k2}(f,g)
Sk_fg = [Tk_fxy Tk_gxy];

% Remove optimal column
Ak_fg = Sk_fg;
Ak_fg(:,idx_col) = [];
ck = Sk_fg(:,idx_col);

% Build the matrix A_{k1,k2}(zf,zg)
Ak_zfzg = zeros(size(Ak_fg));

% Build the vector removed from S_{k}(zf,zg)
hk = zeros(size(ck));

% Initialise vector of perturbations of f and g
z = zeros(nNonZeros_fxy + nNonZeros_gxy,1);

% Get the vector of coefficients zf(x,y)
v_zf = z(1:nNonZeros_fxy );

% Get the vector of coefficeints zg
v_zg = z(nNonZeros_fxy+1:end);

% Get zf as a matrix
mat_zfxy = GetAsMatrix([v_zf; zeros(nZeros_fxy,1)], m, m);

% Get zg as a matrix
mat_zgxy = GetAsMatrix([v_zg; zeros(nZeros_gxy,1)], n, n);

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

Yk = BuildY_TotalDegree_STLN(x,m,n,k);

% Get vector of coefficients of f(x,y)
v_fxy = GetAsVector(fxy_matrix);
v_fxy = v_fxy(1:nNonZeros_fxy);

% Get vector of coefficients of g(x,y)
v_gxy = GetAsVector(gxy_matrix);
v_gxy = v_gxy(1:nNonZeros_gxy);

% % Test the funciton which constructs Y_{k}
test1a = Yk * [v_fxy;v_gxy];
test1b = Ak_fg * xk;
test1 = norm(test1a - test1b);
display(test1);

% Build the matrix P_{t}
% Where P * [f;g] = c_{t}
Pk = BuildP_TotalDegree_STLN(m,n,idx_col,k);

% Test the function which constructs P_{k}
test2a = Pk * [v_fxy;v_gxy];
test2b = ck;
test2 = norm(test2a - test2b);
display(test2)

% % Build the Matrix C consisting of H_{z} and H_{x}
H_z = Yk - Pk;

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
    nNonZeroEntries_z      = nNonZeros_fxy + nNonZeros_gxy;
    delta_zk        = y_lse(1:nNonZeroEntries_z);
    delta_xk        = y_lse((nNonZeroEntries_z+1):end);
    
    % Update z and x
    z       = z + delta_zk;
    xk    = xk + delta_xk;
    
    % Split vector z into vectors z_f and z_g
    v_zf = z(1 : nNonZeros_fxy);
    v_zg = z(nNonZeros_fxy + 1 : end);
    
    % Get coefficients of zf(x,y) in a matrix
    mat_zfxy = GetAsMatrix([v_zf; zeros(nZeros_fxy,1)], m, m);
    
    % Get coefficients of zg(x,y) in a matrix
    mat_zgxy = GetAsMatrix([v_zg;zeros(nZeros_gxy,1)], n, n);
    
    % Build the matrix B_{t} = [E1(zf) E2(zg)]
    E1 = BuildT1_Total(mat_zfxy, m, n-k);
    E2 = BuildT1_Total(mat_zgxy, n, m-k);
    
    % Build the matrix B_{t} equivalent to S_{t}
    Bt = [E1 E2];
    
    % Get the matrix E_{t} with optimal column removed
    Ak_zfzg = Bt;
    Ak_zfzg(:,idx_col) = [];
    
    % Get the column vector h_{t}, the optimal column removed from B_{t},
    % and equivalent to c_{t} removed from S_{t}
    hk = Bt(:,idx_col);
    
    % Get the updated vector x = [x1 x2] where x1 and x2 are vectors.
    % S(x1,x2)*[f;g] = ct
    % x_ls = SolveAx_b(At+Et,ct+ht);
    
    x = [...
        xk(1:idx_col-1);...
        0;...
        xk(idx_col:end)];
    
    % Build the matrix Y_{t} where Y_{t}(x)*z = E_{t}(z) * x
    Yk = BuildY_TotalDegree_STLN(x,m,n,k);
    
    % Get the residual vector
    res_vec = (ck+hk) - ((Ak_fg+Ak_zfzg)*xk);
    
    % Update the matrix C
    H_z = Yk - Pk;
    H_x = Ak_fg + Ak_zfzg;
    C = [H_z H_x];
    
    % Update fnew - used in LSE Problem.
    f = -(yy-start_point);
    
    % Update the termination criterion
    condition(ite) = norm(res_vec)./norm(ck+hk) ;
    
    
end

fprintf([mfilename ' : ' sprintf('Required number of iterations: %i\n',ite)]);

PlotGraphs_STLN();


fxy_lr = fxy_matrix + mat_zfxy;
gxy_lr = gxy_matrix + mat_zgxy;

end


