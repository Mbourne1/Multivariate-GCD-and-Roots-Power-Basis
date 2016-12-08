function [uxy, vxy, wxy] = GetQuotients_Total_3Polys(fxy, gxy, hxy, m, n, o, k)
% Given two polynomials and the knowledge of the degree of the GCD. Obtain
% the two quotient polynomials u and v.
%
% % Inputs.
%
% [fxy, gxy, hxy] : Coefficients of polynomail f(x,y) g(x,y) and h(x,y)
%
% [m, n, o] : Total degree of f(x,y), g(x,y) and h(x,y)
%
% k : Degree of GCD
%
% % Outputs
% 
% [uxy, vxy, wxy] : Coefficients of polynomial u(x,y), v(x,y) and w(x,y)


% %
% %
% Replace fxy with the padded version
% Get degree of f(x,y) with respect to x and y
[m1,m2] = GetDegree(fxy);
padd_matrix = zeros(m+1,m+1);
padd_matrix(1:m1+1,1:m2+1) = fxy;
fxy = padd_matrix;

% %
% %
% Replace gxy with the padded version
% Get degree of g(x,y) with respect to x and y
[n1,n2] = GetDegree(gxy);
padd_matrix = zeros(n+1,n+1);
padd_matrix(1:n1+1,1:n2+1) = gxy;
gxy = padd_matrix;


%
%
[o1, o2] = GetDegree(hxy);
padd_matrix = zeros(o+1,o+1);
padd_matrix(1:o1+1,1:o2+1) = hxy;
hxy = padd_matrix;


% Build the Sylvester Subresultant matrix S
Sk = BuildSylvesterMatrix_Total_3Polys(fxy, gxy, hxy, m, n, o, k);

% Get the optimal column for removal from S(f,g)
opt_col = GetOptimalColumn_Total(Sk);


% % Having found the optimal column, obtain u and v the quotient polynomials.
Akj = Sk;
cki = Sk(:,opt_col);
Akj(:,opt_col) = [];

x_ls = SolveAx_b(Akj,cki);

% Obtain the solution vector x = [-v;u]
vecx =[
    x_ls(1:(opt_col)-1);
    -1;
    x_ls(opt_col:end);
    ]  ;

% Get number of coefficients in u(x,y) and v(x,y)
nCoeffs_vxy = nchoosek(n-k+2,2);
nCoeffs_wxy = nchoosek(o-k+2,2);
nCoeffs_uxy = nchoosek(m-k+2,2);

% get the vector of coefficients of v
vxy_calc = vecx(1:nCoeffs_vxy);

wxy_calc = vecx(nCoeffs_vxy + 1 : nCoeffs_vxy + nCoeffs_wxy);

% get the vector of coefficients of u
uxy_calc = (-1) .* vecx(nCoeffs_vxy + nCoeffs_wxy + 1 : end);
            
        

% %
% %
% Get u(x,y) and v(x,y) in matrix form
zeros_vww = zeros(nchoosek(n-k-1+2,2),1);
zeros_uww = zeros(nchoosek(m-k-1+2,2),1);
zeros_www = zeros(nchoosek(o-k-1+2,2),1);

uxy = GetAsMatrix([uxy_calc;zeros_uww],m-k,m-k);
wxy = GetAsMatrix([wxy_calc;zeros_www],o-k,o-k);
vxy = GetAsMatrix([vxy_calc;zeros_vww],n-k,n-k);





end
