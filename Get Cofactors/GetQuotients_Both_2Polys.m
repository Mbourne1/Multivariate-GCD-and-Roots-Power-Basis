function [uxy_calc_matrix,vxy_calc_matrix] = ...
    GetQuotients_Both_2Polys(fxy, gxy, m, n, k, k1, k2)
% GetQuotients(fxy_matrix,gxy_matrix,t1,t2,opt_alpha,th1,th2)
%
% Given two polynomials and the knowledge of the degree of the GCD. Obtain
% the two quotient polynomials u and v.
%
% % Inputs.
%
% fxy : Coefficients of the polynomial f(x,y).
%
% gxy : coefficients of the polynomial g(x,y).
%
% k : Total degree of d(x,y)
%
% k1 : Degree of GCD d(x,y).
%
% k2 : Degree of GCD d(x,y).

% Get the degree of polynomial f(x,y).
[m1,m2] = GetDegree(fxy);

% Get the degree of polynomial g(x,y).
[n1,n2] = GetDegree(gxy);

% Build the Sylvester matrix S_{k,k1,k2}
Skk1k2 = BuildSylvesterMatrix_Both_2Polys(fxy,gxy,m,n,k,k1,k2);

% % Get the optimal column for removal.
idx_col = GetOptimalColumn_Both(Skk1k2);

% Having found the optimal column, obtain u and v the quotient polynomials.
Atj = Skk1k2;
cki = Skk1k2(:,idx_col);
Atj(:,idx_col) = [];

% Get the solution vector.
x_ls = SolveAx_b(Atj,cki);

% Get number of non-zero entries in u(x,y) and v(x,y)
nNoneZeros_uxy = GetNumNonZeros(m1-k1,m2-k2,m-k);
nNoneZeros_vxy = GetNumNonZeros(n1-k1,n2-k2,n-k);

% Get number of zero entries in u(x,y) and v(x,y)
nZeros_uxy = (m1-k1+1) * (m2-k2+1) - nNoneZeros_uxy;
nZeros_vxy = (n1-k1+1) * (n2-k2+1) - nNoneZeros_vxy;

% Obtain the solution vector x = [-v;u]
vecx =[
    x_ls(1:(idx_col)-1);
    -1;
    x_ls(idx_col:end);
    ]  ;

% Get the vector of coefficients of v(x,y)
vxy_calc = [...
            vecx(1:(nNoneZeros_vxy));
            zeros(nZeros_vxy,1)
          ];
      
% Get the vector of coefficients of u(x,y)
uxy_calc = [...
            (-1).*vecx( nNoneZeros_vxy + 1 :end);
            zeros(nZeros_uxy,1);
            ];
        



% Arrange u(x,y) as a matrix.
uxy_calc_matrix = GetAsMatrix(uxy_calc,m1-k1,m2-k2);

% Arrange v(x,y) as a matrix.
vxy_calc_matrix = GetAsMatrix(vxy_calc,n1-k1,n2-k2);





end
