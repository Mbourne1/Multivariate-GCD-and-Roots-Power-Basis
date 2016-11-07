function [fxy_lr,gxy_lr] = STLN_Total(fxy_matrix,gxy_matrix,m,n,k,idx_col)
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
nCoeff_f = nchoosek(m+2,2);
nZeros_f = nchoosek(m+1,2);

% Get the number of coefficients in the polynomial g(x,y)
nCoeff_g = nchoosek(n+2,2);
nZeros_g = nchoosek(m+1,2);

% Get the number of coefficients in v(x,y)
nCoeff_v = nchoosek(n-k+2,2);

% Get the number of coefficients in u(x,y)
nCoeff_u = nchoosek(m-k+2,2);

% Build the Sylvester Matrix S_{k1,k2}(f,g)

% Build the partiton T_{n1-k1,n2-k2}(f)
Tk_f = BuildT1_Total(fxy_matrix,m,n-k);

% Build the partition T_{m1-k1,m2-k2}(g)
Tk_g = BuildT1_Total(gxy_matrix,n,m-k);

% Build the Sylvester subresultant matrix S_{k1,k2}(f,g)
Sk_fg = [Tk_f Tk_g];

% Remove optimal column
Ak_fg = Sk_fg;
Ak_fg(:,idx_col) = [];
ct = Sk_fg(:,idx_col);

% Build the matrix A_{k1,k2}(zf,zg)
Ak_zfzg = zeros(size(Ak_fg));

% Build the vector removed from S_{k}(zf,zg)
ht = zeros(size(ct));

% Initialise vector of perturbations of f and g
z = zeros(nCoeff_f + nCoeff_g,1);

% Get the vector of coefficients zf(x,y)
v_zf = z(1:nCoeff_f );

% Get the vector of coefficeints zg
v_zg = z(nCoeff_f+1:end);

% Get zf as a matrix
% EDIT 10/03/2016 - vZ_fxy has zeros removed, so include the zeros to form
% a matrix m1+1 * m2+1

mat_sfxy = GetAsMatrix([v_zf; zeros(nZeros_f,1)], m, m);

% Get zg as a matrix
mat_zgxy = GetAsMatrix([v_zg; zeros(nZeros_g,1)], n, n);

% Get the vector x
% A_{t} x = c_{t}
x_ls = SolveAx_b(Ak_fg,ct);


x = ...
    [
    x_ls(1:idx_col-1);
    0;
    x_ls(idx_col:end);
    ];



% Build the matrix Y_{t}
% Where Y(x) * z = E(z) * x
% #########################################################################
% To Do - Write this function

Yk = BuildY_TotalDegree_STLN(x,m,n,k);
% #########################################################################


% Get vector of coefficients of f(x,y)
%%v_fxy = GetAsVector(fxy_matrix);

% Get vector of coefficients of g(x,y)
%%v_gxy = GetAsVector(gxy_matrix);

% Test
%%test1 = Yk * [v_fxy;v_gxy];
%%test2 = Ak_fg * x_ls;
%%norm(test1-test2)

% Build the matrix P_{t}
% Where P * [f;g] = c_{t}
Pk = BuildP_TotalDegree_STLN(m,n,idx_col,k);

% Test
%%test1 = P * [v_fxy;v_gxy];
%%test2 = ct;
%%norm(test1-test2)


% Get initial residual (A_{t}+E_{t})x = (c_{t} + h_{t})
x_ls = SolveAx_b(Ak_fg+Ak_zfzg,ct+ht);

% % Build the Matrix C consisting of H_{z} and H_{x}
H_z = Yk - Pk;

H_x = Ak_fg + Ak_zfzg;

C = [H_z H_x];

E = blkdiag( eye(nCoeff_f + nCoeff_g) , eye(nCoeff_u + nCoeff_v - 1));


% Get initial residual (A_{t}+E_{t})x = (c_{t} + h_{t})
res_vec = (ct + ht) - (Ak_fg*x_ls);

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
    nEntries_z      = nCoeff_f + nCoeff_g;
    delta_zk        = y_lse(1:nEntries_z);
    delta_xk        = y_lse((nEntries_z+1):end);
    
    % Update z and x
    z       = z + delta_zk;
    x_ls    = x_ls + delta_xk;
    
    % Split vector z into vectors z_f and z_g
    v_zf = z(1 : nCoeff_f);
    v_zg = z(nCoeff_f + 1 : end);
    
    % Get zf as a matrix
    mat_sfxy = GetAsMatrix([v_zf; zeros(nCoeff_f,1)], m, m);
    
    % Get zg as a matrix
    mat_zgxy = GetAsMatrix([v_zg;zeros(nCoeff_g,1)], n, n);
    
    % Build the matrix B_{t} = [E1(zf) E2(zg)]
    E1 = BuildT1_Total(mat_sfxy, m, n-k);
    E2 = BuildT1_Total(mat_zgxy, n, m-k);
    
    % Build the matrix B_{t} equivalent to S_{t}
    Bt = [E1 E2];
    
    % Get the matrix E_{t} with optimal column removed
    Ak_zfzg = Bt;
    Ak_zfzg(:,idx_col) = [];
    
    % Get the column vector h_{t}, the optimal column removed from B_{t},
    % and equivalent to c_{t} removed from S_{t}
    ht = Bt(:,idx_col);
    
    % Get the updated vector x = [x1 x2] where x1 and x2 are vectors.
    % S(x1,x2)*[f;g] = ct
    % x_ls = SolveAx_b(At+Et,ct+ht);
    
    x = [...
        x_ls(1:idx_col-1);...
        0;...
        x_ls(idx_col:end)];
    
    % Build the matrix Y_{t} where Y_{t}(x)*z = E_{t}(z) * x
    Yk = BuildY_TotalDegree_STLN(x,m,n,k);
    
    % Get the residual vector
    res_vec = (ct+ht) - ((Ak_fg+Ak_zfzg)*x_ls);
    
    % Update the matrix C
    H_z = Yk - Pk;
    H_x = Ak_fg + Ak_zfzg;
    C = [H_z H_x];
    
    % Update fnew - used in LSE Problem.
    f = -(yy-start_point);
    
    % Update the termination criterion
    condition(ite) = norm(res_vec)./norm(ct+ht) ;
    
    
end

fprintf([mfilename ' : ' sprintf('Required number of iterations: %i\n',ite)]);

PlotGraphs_STLN();

fxy_lr = fxy_matrix + mat_sfxy;
gxy_lr = gxy_matrix + mat_zgxy;

display([GetAsVector(fxy_matrix) GetAsVector(mat_sfxy)]);
display([GetAsVector(gxy_matrix) GetAsVector(mat_zgxy)]);

end


