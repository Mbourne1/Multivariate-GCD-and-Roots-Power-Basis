function [uxy, vxy, wxy] = GetQuotients_Relative_3Polys(fxy, gxy, hxy, k1, k2)
% GetQuotients(fxy,gxy,k1,k2)
%
% Given two polynomials and the knowledge of the degree of the GCD. Obtain
% the two quotient polynomials u and v.
%
% % Inputs.
%
% [fxy, gxy, hxy] : Coefficients of the polynomial f(x,y), g(x,y) and
% h(x,y)
%
% t1 : Degree of GCD d(x,y) with respect to x
%
% t2 : Degree of GCD d(x,y) with respect to y
%
% % Outputs
%
% [uxy, vxy, wxy]
%

% Get the degree of polynomial f(x,y).
[m1, m2] = GetDegree_Bivariate(fxy);

% Get the degree of polynomial g(x,y).
[n1, n2] = GetDegree_Bivariate(gxy);

% Get the degree of polynomial h(x,y)
[o1, o2] = GetDegree_Bivariate(hxy);

Sk = BuildSylvesterMatrix_Relative_3Polys(fxy, gxy, hxy, k1, k2);

% % Get the optimal column for removal
opt_col = GetOptimalColumn_Relative(Sk);

% Having found the optimal column, obtain u and v the quotient polynomials.
Atj = Sk;
cki = Sk(:,opt_col);
Atj(:,opt_col) = [];

% Get the solution vector.
x_ls = SolveAx_b(Atj,cki);

% Obtain the solution vector x = [-v;u]
vecx =[
    x_ls(1:(opt_col)-1);
    -1;
    x_ls(opt_col:end);
    ]  ;

% Get number of coefficients in u(x,y) and v(x,y) and w(x,y)
nCoeffs_vxy = (n1-k1+1) * (n2-k2+1);
nCoeffs_wxy = (o1-k1+1) * (o2-k2+1);

% Get the vector of coefficients of v
vxy_calc = vecx(1:nCoeffs_vxy);

wxy_calc = vecx(nCoeffs_vxy + 1 : nCoeffs_vxy + nCoeffs_wxy);
      
% Get the vector of coefficients of u
uxy_calc = (-1).*vecx(nCoeffs_vxy + nCoeffs_wxy +1 : end);
        

% % Get u and v in matrix form
% Arrange uw into a matrix form based on their dimensions.
uxy = GetAsMatrix(uxy_calc, m1-k1, m2-k2);
wxy = GetAsMatrix(wxy_calc, o1-k1, o2-k2);
vxy = GetAsMatrix(vxy_calc, n1-k1, n2-k2);


end
