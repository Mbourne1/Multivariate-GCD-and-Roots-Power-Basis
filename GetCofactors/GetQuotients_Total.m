function [uxy_calc_mtrx,vxy_calc_mtrx] = ...
    GetQuotients_Total(fxy_matrix,gxy_matrix,m,n,t)
% Given two polynomials and the knowledge of the degree of the GCD. Obtain
% the two quotient polynomials u and v.
%
% % Inputs.
%
% fxy_matrix : Coefficients of polynomail f(x,y)
%
% gxy_matrix : Coefficinets of polynomial g(x,y)
%
% m : Total degree of f(x,y)
%
% n : Total degree of g(x,y)
%
% t : Degree of GCD


% %
% %
% Replace fxy_matrix with the padded version
% Get degree of f(x,y) with respect to x and y
[m1,m2] = GetDegree(fxy_matrix);
padd_matrix = zeros(m+1,m+1);
padd_matrix(1:m1+1,1:m2+1) = fxy_matrix;
fxy_matrix = padd_matrix;

% %
% %
% Replace gxy_matrix with the padded version
% Get degree of g(x,y) with respect to x and y
[n1,n2] = GetDegree(gxy_matrix);
padd_matrix = zeros(n+1,n+1);
padd_matrix(1:n1+1,1:n2+1) = gxy_matrix;
gxy_matrix = padd_matrix;

% % 
% %
% Build the partitions of the Sylvester matrix S_{t}

% Build the first partition containing coefficients of fxy
T1 = BuildT1_TotalDegree(fxy_matrix,m,n-t);

% Build the second partition containing coefficients of gxy
T2 = BuildT1_TotalDegree(gxy_matrix,n,m-t);

% Concatenate the two partitions
St = [T1 T2];

% %
% % 
% Get the optimal column for removal from S(f,g)
opt_col = GetOptimalColumn_Total(fxy_matrix,gxy_matrix,m,n,t);


% % Having found the optimal column, obtain u and v the quotient polynomials.
Atj = St;
cki = St(:,opt_col);
Atj(:,opt_col) = [];

x_ls = SolveAx_b(Atj,cki);

% Obtain the solution vector x = [-v;u]
vecx =[
    x_ls(1:(opt_col)-1);
    -1;
    x_ls(opt_col:end);
    ]  ;

% Get number of coefficients in u(x,y) and v(x,y)
num_coeff_v = nchoosek(n-t+2,2);
num_coeff_u = nchoosek(m-t+2,2);

% get the vector of coefficients of v
vxy_calc = vecx(1:num_coeff_v);
      
% get the vector of coefficients of u
uxy_calc = (-1).*vecx(num_coeff_v+1:end);
            
        

% %
% %
% Get u(x,y) and v(x,y) in matrix form
zeros_vww = zeros(nchoosek(n-t-1+2,2),1);
zeros_uww = zeros(nchoosek(m-t-1+2,2),1);

uxy_calc_mtrx = GetAsMatrix([uxy_calc;zeros_uww],m-t,m-t);
vxy_calc_mtrx = GetAsMatrix([vxy_calc;zeros_vww],n-t,n-t);





end
