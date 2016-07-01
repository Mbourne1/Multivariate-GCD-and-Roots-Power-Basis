function [uxy_calc_matrix,vxy_calc_matrix] = ...
    GetQuotients_both(fxy_matrix,gxy_matrix,m,n,t,t1,t2,opt_alpha,th1,th2)
% GetQuotients(fxy_matrix,gxy_matrix,t1,t2,opt_alpha,th1,th2)
%
% Given two polynomials and the knowledge of the degree of the GCD. Obtain
% the two quotient polynomials u and v.
%
%       Inputs.
%
%   fxy_matrix : Coefficients of the polynomial f(x,y).
%
%   gxy_matrix : coefficients of the polynomial g(x,y).
%
%   t1 : Degree of GCD d(x,y).
%
%   t2 : Degree of GCD d(x,y).
%
%   opt_alpha : alpha.
%
%   th1 : \theta_{1}.
%   
%   th2 : \theta_{2}.


% Get the degree of polynomial f(x,y).
[m1,m2] = GetDegree(fxy_matrix);

% Get the degree of polynomial g(x,y).
[n1,n2] = GetDegree(gxy_matrix);

% % Preprocess

% Get f(w,w)
fww_matrix = GetWithThetas(fxy_matrix,th1,th2);

% Get g(w,w) 
gww_matrix = GetWithThetas(gxy_matrix,th1,th2);

% % Build the partitions of the Sylvester matrix S_{t}

% Build the first partition containing coefficients of fxy
T1 = BuildT1(fww_matrix,n1-t1,n2-t2);

% Build the second partition containing coefficients of gxy
T2 = BuildT1(opt_alpha.*gww_matrix,m1-t1,m2-t2);

% Get number of non-zero entries in u(x,y) and v(x,y)
nNoneZeros_uxy = GetNumNonZeros(m1-t1,m2-t2,m-t);
nNoneZeros_vxy = GetNumNonZeros(n1-t1,n2-t2,n-t);

% Get number of zero entries in u(x,y) and v(x,y)
nZeros_uxy = (m1-t1+1) * (m2-t2+1) - nNoneZeros_uxy;
nZeros_vxy = (n1-t1+1) * (n2-t2+1) - nNoneZeros_vxy;

% %
% %
% Remove Columns from S

% Remove the zero columns from T_{n1-t1,n2-t2}
T1 = T1(:,1:nNoneZeros_vxy);

% Remove the zero columns from T_{m1-t1,m2-t2}
T2 = T2(:,1:nNoneZeros_uxy);

% %
% % 
% Remove Rows from S

nNoneZeros_fv = GetNumNonZeros(m1+n1-t1,m2+n2-t2,m+n-t);
nNoneZeros_gu = GetNumNonZeros(n1+m1-t1,n2+m2-t2,n+m-t);

T1 = T1(1:nNoneZeros_fv,:);
T2 = T2(1:nNoneZeros_gu,:);


% Form the Sylvester matrix.
St = [T1 T2];

% % Get the optimal column for removal.
opt_col = GetOptimalColumn_both(fww_matrix,opt_alpha.*gww_matrix,m,n,t,t1,t2);


% Having found the optimal column, obtain u and v the quotient polynomials.
Atj = St;
cki = St(:,opt_col);
Atj(:,opt_col) = [];

% Get the solution vector.
x_ls = SolveAx_b(Atj,cki);

% Obtain the solution vector x = [-v;u]
vecx =[
    x_ls(1:(opt_col)-1);
    -1;
    x_ls(opt_col:end);
    ]  ;

vecx = vecx./vecx(1);

% Get the vector of coefficients of v
vww_calc = [...
            vecx(1:(nNoneZeros_vxy));
            zeros(nZeros_vxy,1)
          ];
      
% Get the vector of coefficients of u
uww_calc = [...
            (-1).*vecx( nNoneZeros_vxy + 1 :end);
            zeros(nZeros_uxy,1);
            ];
        



% Arrange u(w,w) as a matrix.
uww_calc_matrix = GetAsMatrix(uww_calc,m1-t1,m2-t2);

% Arrange v(w,w) as a matrix.
vww_calc_matrix = GetAsMatrix(vww_calc,n1-t1,n2-t2);

% % Get u(x,y) and v(x,y) from u(w,w) and v(w,w)
uxy_calc_matrix = GetWithoutThetas(uww_calc_matrix,th1,th2);
vxy_calc_matrix = GetWithoutThetas(vww_calc_matrix,th1,th2);






end
