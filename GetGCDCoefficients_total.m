function dxy_matrix = GetGCDCoefficients_total(fxy,gxy,...
    uxy, vxy,alpha,th1,th2,m,n,t)
% Given the matrices of coefficients of f(x,y) and g(x,y), the quotient
% polynomials u(x,y) and v(x,y), and optimal values for alpha, theta_{1}
% and theta_{2}, calculate the coefficients of the GCD d(x,y).
%
% %         Inputs
%
% fxy : Coefficients of polynomial f(x,y)
%
% gxy : Coefficients of polynomial g(x,y)
%
% uxy : Coefficients of polynomial u(x,y)
%
% vxy : Coefficients of polynomial v(x,y)
%
% alpha : Optimal value of alpha
%
% th1 : Optimal value of theta_{1}
%
% th2 : Optimal value of theta_{2}


% % Preprocess u and v

% Preprocess u(x,y) to obtain u(w,w)
uww_matrix = GetWithThetas(uxy,th1,th2);

% Preprocess v(x,y) to obtain v(w,w)
vww_matrix = GetWithThetas(vxy,th1,th2);

% % Build Matrix C

% Build the Cauchy matrix of coefficients of u(w,w)
C1 = BuildT1_TotalDegree(uww_matrix,m,n-t);

% Build the Cauchy matrix of coefficients of v(w,w)
C2 = BuildT1_TotalDegree(vww_matrix,n,m-t);

% Build the RHS vector of coefficients of f and g
C = [C1;C2];

%% Build Vector f(w,w)
% Get f with thetas included

% padd f so that it is in terms of its total degree

[m1,m2] = GetDegree(fxy);
padd = zeros(m+1,m+1);
padd(1:m1+1,1:m2+1) = fxy;
fxy = padd;

[n1,n2] = GetDegree(gxy);
padd = zeros(n+1,n+1);
padd(1:n1+1,1:n2+1) = gxy;
gxy = padd;

% Get f(w,w) from f(x,y)
fww_matrix = GetWithThetas(fxy,th1,th2);

% Get fww_matrix as a vector
fww_vec = GetAsVector(fww_matrix);

% Remove the zeros associated with the polynomial by total degree
num_elements_f = nchoosek(m+2,2);
fww_vec = fww_vec(1:num_elements_f);

% % Build Vector g(w,w)
gww_matrix = GetWithThetas(gxy,th1,th2);

% get gww_matrix as a vector
gww_vec = GetAsVector(gww_matrix);

% Remove the zeros associated with the polynomial by total degree
num_elements_g = nchoosek(n+2,2);
gww_vec = gww_vec(1:num_elements_g);

% % Build the RHS vector
rhs_vec = [fww_vec;
           alpha .*gww_vec];

% % Calculate the solution vector

% Calculate the x vector by pinv 
x = SolveAx_b(C,rhs_vec);
    
% Get the residual associated with the solution x. (Small residual implies good approximation)    
residual = pinv(C)*rhs_vec - x;

% Padd d(w,w) with zeros so that it can be put back into matrix form
dww_vec = [dww_vec ; zeros(nchoosek(t+2-1,2),1)];

% Arrange dw into a matrix form based on its dimensions.
dww_calc_matrix = getAsMatrix(dww_vec,t,t);

% % Obtain d(x,y) from d(w,w)
dxy_matrix = GetWithoutThetas(dww_calc_matrix,th1,th2);

end