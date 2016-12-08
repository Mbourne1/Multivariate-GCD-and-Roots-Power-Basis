function [uxy, vxy] = GetQuotients_Relative_2Polys(fxy, gxy, k1, k2)
% GetQuotients(fxy,gxy,k1,k2)
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
% t1 : Degree of GCD d(x,y).
%
% t2 : Degree of GCD d(x,y).


% Get the degree of polynomial f(x,y).
[m1,m2] = GetDegree(fxy);

% Get the degree of polynomial g(x,y).
[n1,n2] = GetDegree(gxy);

Sk = BuildSylvesterMatrix_Relative_2Polys(fxy, gxy, k1, k2);

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

% Get number of coefficients in u(x,y) and v(x,y)
num_coeff_v = (n1-k1+1) * (n2-k2+1);

% Get the vector of coefficients of v
vxy_calc = vecx(1:num_coeff_v);
      
% Get the vector of coefficients of u
uxy_calc = (-1).*vecx(num_coeff_v+1:end);
        

% % Get u and v in matrix form
% Arrange uw into a matrix form based on their dimensions.
uxy = GetAsMatrix(uxy_calc, m1-k1, m2-k2);
vxy = GetAsMatrix(vxy_calc, n1-k1, n2-k2);


end
