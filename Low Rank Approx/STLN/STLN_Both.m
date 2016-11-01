function [fxy_lr,gxy_lr] = STLN_Both(fxy_matrix,gxy_matrix,m,n,t,t1,t2,idx_opt_col)
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
Tf = BuildT1(fxy_matrix,n1-t1,n2-t2);
Tg = BuildT1(gxy_matrix,m1-t1,m2-t2);

% %
% %
% Remove the columns of T(f) and T(g) which correspond to the zeros in u(x,y)
% and v(x,y) which are removed from the solution vector x.
Tf = Tf(:,1:nNonZeros_vxy);
Tg = Tg(:,1:nNonZeros_uxy);

% % 
% %
% Remove the extra rows of T1 and T2 associated with zeros of f*v and g*u
Tf = Tf(1:nNonZeros_fv,:);
Tg = Tg(1:nNonZeros_gu,:);

% Build the matrix S_{t_{1},t_{2}} with reduced columns
St = [Tf Tg];

% Remove optimal column
At = St;
At(:,idx_opt_col) = [];
ct = St(:,idx_opt_col);

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

x = ...
    [
    x_ls(1:idx_opt_col-1);
    0;
    x_ls(idx_opt_col:end);
    ];


% Build the matrix Y_{t}
% Where Y(x) * z = E(z) * x
Yt = BuildY_BothDegree_STLN(x,m,m1,m2,n,n1,n2,t,t1,t2);


v_fxy = GetAsVector(fxy_matrix);
v_fxy = v_fxy(1:nNonZeros_fxy,:);

v_gxy = GetAsVector(gxy_matrix);
v_gxy = v_gxy(1:nNonZeros_gxy,:);


test1 = Yt * [v_fxy;v_gxy];
test2 = At * x_ls;
norm(test1-test2)

% Build the matrix P_{t}
% Where P * [f;g] = c_{t}
Pt = BuildP_BothDegree_STLN(m,m1,m2,n,n1,n2,t,t1,t2,idx_opt_col);
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
    Et(:,idx_opt_col) = [];
    
    % Get the column vector h_{t}, the optimal column removed from B_{t},
    % and equivalent to c_{t} removed from S_{t}
    ht = Bt(:,idx_opt_col);
    
    % Get the updated vector x
    %x_ls = SolveAx_b(At + Et,ct + ht);
    
    x = [...
        x_ls(1:idx_opt_col-1);...
        0;...
        x_ls(idx_opt_col:end)];
    
    % Build the matrix Y_{t} where Y_{t}(x)*z = E_{t}(z) * x
    Yt = BuildY_BothDegree_STLN(x,m,m1,m2,n,n1,n2,t,t1,t2);
    
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


%if condition(ite) < condition(1)
    
    fxy_lr = fxy_matrix + matZ_fxy;
    gxy_lr = gxy_matrix + matZ_gxy;
%else
%    fxy_lr = fxy_matrix;
%    gxy_lr = gxy_matrix;
%end

PlotGraphs_STLN


end

