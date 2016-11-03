function [uxy,vxy] = GetQuotients_Total(fxy,gxy,m,n,k)
% Given two polynomials and the knowledge of the degree of the GCD. Obtain
% the two quotient polynomials u and v.
%
% % Inputs.
%
% fxy: Coefficients of polynomail f(x,y)
%
% gxy: Coefficinets of polynomial g(x,y)
%
% m : Total degree of f(x,y)
%
% n : Total degree of g(x,y)
%
% k : Degree of GCD
%
% % Outputs
% 
% uxy : Coefficients of polynomial u(x,y)
%
% vxy : Coefficeints of polynomial v(x,y)


% %
% %
% Replace fxy_matrix with the padded version
% Get degree of f(x,y) with respect to x and y
[m1,m2] = GetDegree(fxy);
padd_matrix = zeros(m+1,m+1);
padd_matrix(1:m1+1,1:m2+1) = fxy;
fxy = padd_matrix;

% %
% %
% Replace gxy_matrix with the padded version
% Get degree of g(x,y) with respect to x and y
[n1,n2] = GetDegree(gxy);
padd_matrix = zeros(n+1,n+1);
padd_matrix(1:n1+1,1:n2+1) = gxy;
gxy = padd_matrix;


% Build the Sylvester Subresultant matrix S
Sk = BuildSylvesterMatrix_Total(fxy,gxy,m,n,k);

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
num_coeff_v = nchoosek(n-k+2,2);
num_coeff_u = nchoosek(m-k+2,2);

% get the vector of coefficients of v
vxy_calc = vecx(1:num_coeff_v);
      
% get the vector of coefficients of u
uxy_calc = (-1).*vecx(num_coeff_v+1:end);
            
        

% %
% %
% Get u(x,y) and v(x,y) in matrix form
zeros_vww = zeros(nchoosek(n-k-1+2,2),1);
zeros_uww = zeros(nchoosek(m-k-1+2,2),1);

uxy = GetAsMatrix([uxy_calc;zeros_uww],m-k,m-k);
vxy = GetAsMatrix([vxy_calc;zeros_vww],n-k,n-k);





end
