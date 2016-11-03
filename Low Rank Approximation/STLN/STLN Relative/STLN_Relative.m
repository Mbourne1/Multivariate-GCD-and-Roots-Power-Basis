function [fxy_lr,gxy_lr] = STLN_Relative(fxy_matrix,gxy_matrix,t1,t2,idx_col)
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
% t1 : Degree of d(x,y) with respect to x
%
% t2 : Degree of d(x,y) with respect to y
%
% idx_col : Index of optimal column for removal from S_{t_{1},t_{2}}(f,g)
%
% Outputs 
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
nCoeff_f = (m1+1) * (m2+1);

% Get the number of coefficients in the polynomial g(x,y)
nCoeff_g = (n1+1) * (n2+1);

% Get the number of coefficients in v(x,y)
nCoeff_v = (n1-t1+1) * (n2-t2+1);

% Get the number of coefficients in u(x,y)
nCoeff_u = (m1-t1+1) * (m2-t2+1);

% Build the Sylvester Matrix S(f,g)
Tf = BuildT1_Relative(fxy_matrix,n1-t1,n2-t2);
Tg = BuildT1_Relative(gxy_matrix,m1-t1,m2-t2);
St = [Tf Tg];

% Remove optimal column
At = St;
At(:,idx_col) = [];
ct = St(:,idx_col);

% Build the matrix Et
Et = zeros(size(At));
ht = zeros(size(ct));

% EDIT 10/03/2016
z = zeros(nCoeff_f + nCoeff_g,1);

% Get the vector of coefficients zf
vZ_fxy = z(1:nCoeff_f );

% Get the vector of coefficeints zg
vZ_gxy = z(nCoeff_f+1:end);

% Get zf as a matrix
% EDIT 10/03/2016 - vZ_fxy has zeros removed, so include the zeros to form
% a matrix m1+1 * m2+1
matZ_fxy = GetAsMatrix(vZ_fxy, m1, m2);

% Get zg as a matrix
matZ_gxy = GetAsMatrix(vZ_gxy, n1, n2);

% Get the vector x
% A_{t} x = c_{t}
x_ls = SolveAx_b(At,ct);

display(x_ls);

x = ...
    [
    x_ls(1:idx_col-1);
    0;
    x_ls(idx_col:end);
    ];



% Build the matrix Y_{t}
% Where Y(x) * z = E(z) * x
Y = BuildY_RelativeDegree_STLN(x,m1,m2,n1,n2,t1,t2);

v_fxy = GetAsVector(fxy_matrix);

v_gxy = GetAsVector(gxy_matrix);

display(Y)
test1 = Y * [v_fxy;v_gxy];
test2 = At * x_ls;
norm(test1-test2)

% Build the matrix P_{t}
% Where P * [f;g] = c_{t}
P = BuildP_RelativeDegree_STLN(m1,m2,n1,n2,idx_col,t1,t2);
test1 = P * [v_fxy;v_gxy];
test2 = ct;
norm(test1-test2)


% Get initial residual (A_{t}+E_{t})x = (c_{t} + h_{t})
x_ls = SolveAx_b(At+Et,ct+ht);

% % Build the Matrix C consisting of H_{z} and H_{x}
H_z = Y - P;

H_x = At + Et;

C = [H_z H_x];

E = blkdiag( eye(nCoeff_f + nCoeff_g) , zeros(nCoeff_u + nCoeff_v - 1));


% Get initial residual (A_{t}+E_{t})x = (c_{t} + h_{t})
res_vec = (ct + ht) - At*x_ls;

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
    vZ_fxy = z(1 : nCoeff_f);
    vZ_gxy = z(nCoeff_f + 1 : end);
    
    % Get zf as a matrix
    matZ_fxy = GetAsMatrix(vZ_fxy, m1, m2);
    
    % Get zg as a matrix
    matZ_gxy = GetAsMatrix(vZ_gxy, n1, n2);
    
    % Build the matrix B_{t} = [E1(zf) E2(zg)]
    E1 = BuildT1(matZ_fxy,n1-t1,n2-t2);
    E2 = BuildT1(matZ_gxy,m1-t1,m2-t2);
    
    % Build the matrix B_{t} equivalent to S_{t}
    Bt = [E1 E2];
    
    % Get the matrix E_{t} with optimal column removed
    Et = Bt;
    Et(:,idx_col) = [];
    
    % Get the column vector h_{t}, the optimal column removed from B_{t},
    % and equivalent to c_{t} removed from S_{t}
    ht = Bt(:,idx_col);
    
    % Get the updated vector x
    % x_ls = SolveAx_b(At+Et,ct+ht);
    
    x = [...
        x_ls(1:idx_col-1);...
        0;...
        x_ls(idx_col:end)];
    
    % Build the matrix Y_{t} where Y_{t}(x)*z = E_{t}(z) * x
    Y = BuildY_RelativeDegree_STLN(x,m1,m2,n1,n2,t1,t2);
       
    % Get the residual vector
    res_vec = (ct+ht) - ((At+Et)*x_ls);
    
    % Update the matrix C
    H_z = Y - P;
    H_x = At + Et;
    C = [H_z H_x];
    
    % Update fnew - used in LSE Problem.
    f = -(yy-start_point);
    
    % Update the termination criterion
    condition(ite) = norm(res_vec)./norm(ct+ht) ;
    
    
end

fprintf([mfilename ' : ' sprintf('Required number of iterations: %i\n',ite)]);

PlotGraphs_STLN();

%if condition(ite) < condition(1)
    
    fxy_lr = fxy_matrix + matZ_fxy;
    gxy_lr = gxy_matrix + matZ_gxy;
%else
%    fxy_lr = fxy_matrix;
%    gxy_lr = gxy_matrix;
%end

display([GetAsVector(fxy_matrix) GetAsVector(matZ_fxy)]);
display([GetAsVector(gxy_matrix) GetAsVector(matZ_gxy)]);

end


