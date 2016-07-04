function dxy_matrix = GetGCDCoefficients_both(fxy_matrix,gxy_matrix,...
    uxy_matrix, vxy_matrix,m,n,t,...
    opt_alpha,th1,th2)
% Given the matrices of coefficients of f(x,y) and g(x,y), the quotient
% polynomials u(x,y) and v(x,y), and optimal values for alpha, theta_{1}
% and theta_{2}, calculate the coefficients of the GCD d(x,y).
%
% %         Inputs
%
% fxy_matrix : Coefficients of polynomial f(x,y)
%
% gxy_matrix : Coefficients of polynomial g(x,y)
%
% uxy_matrix : Coefficients of polynomial u(x,y)
%
% vxy_matrix : Coefficients of polynomial v(x,y)
%
% m : Total degree of f(x,y)
% 
% n : Total degree of g(x,y)
% 
% t : Total degree of d(x,y)
%
% opt_alpha : Optimal value alpha
%
% th1 : \theta_{2} : Optimal value theta_{1}
%
% th2 : \theta_{1} : Optimal value theta_{2}

% Calculate the GCD of two bivariate polynomials f(x,y) and g(x,y)

% %
% %
% Get Degree structures

% Get degrees of polynomial f(x,y)
[m1,m2] = GetDegree(fxy_matrix);

% Get degree of polynomial g(x,y)
[n1,n2] = GetDegree(gxy_matrix);

% Get degrees of polynomial u(x,y)
[m1_t1,m2_t2] = GetDegree(uxy_matrix);

% Get degrees of polynomial d(x,y)
t1 = m1-(m1_t1);
t2 = m2-(m2_t2);

% %
% %
% Preprocess u and v

% Preprocess u(x,y) to obtain u(w,w)
uww_matrix = GetWithThetas(uxy_matrix,th1,th2);

% Preprocess v(x,y) to obtain v(w,w)
vww_matrix = GetWithThetas(vxy_matrix,th1,th2);

% %
% %
% Build matrices C(u) and C(v)

% Build the Cauchy matrix of coefficients of u(w,w)
C1 = BuildT1(uww_matrix,t1,t2);

% Build the Cauchy matrix of coefficients of v(w,w)
C2 = BuildT1(vww_matrix,t1,t2);

% % 
% % 
% Remove the columns of C(u) and C(v) corresponding to zeros in vector of 
% coefficients of d(x,y)

% Get number of zeros in d(x,y)
nNonZeros_dxy = GetNumNonZeros(t1,t2,t);
nZeros_dxy = (t1+1) * (t2+1) - nNonZeros_dxy;
C1 = C1(:,1:nNonZeros_dxy);
C2 = C2(:,1:nNonZeros_dxy);

% % 
% %
% Remove the rows of C(u) and C(v) corresponding to the zeros in
% u(x,y)*d(x,y) and v(x,y)*d(x,y) which are the zeros in f(x,y) and g(x,y)
% respectively.
nNonZeros_ud = GetNumNonZeros(m1,m2,m);
nNonZeros_vd = GetNumNonZeros(n1,n2,n);
C1 = C1(1:nNonZeros_ud,:);
C2 = C2(1:nNonZeros_vd,:);

C = [C1; C2];

% % 
% % 
% Preprocess f(x,y) and g(x,y) and get in vector form

% Get f(x,y) with thetas included f(w,w)
fww_matrix = GetWithThetas(fxy_matrix,th1,th2);

% Get fww_matrix as a vector
fww_vec = GetAsVector(fww_matrix);

% Remove the zeros from f(x,y) 
fww_vec = fww_vec(1:nNonZeros_ud,:);

% Get g(x,y) with thetas included g(w,w)
gww_matrix = GetWithThetas(gxy_matrix,th1,th2);

% get gww_matrix as a vector
gww_vec = GetAsVector(gww_matrix);

% Remove the zeros from g(x,y)
gww_vec = gww_vec(1:nNonZeros_vd,:);


% Build the RHS vector
rhs_vec = [fww_vec;
           opt_alpha .*gww_vec];

%% Calculate the solution vector

% Calculate the x vector by pinv       
x = SolveAx_b(C,rhs_vec);


% Append the removed zeros
dww_vec = ...
    [
        x;
        zeros(nZeros_dxy,1);
        ];

% Arrange dw into a matrix form based on its dimensions.
dww_calc_mtrx = GetAsMatrix(dww_vec,t1,t2);

% % Obtain d(x,y) from d(w,w)
dxy_matrix = GetWithoutThetas(dww_calc_mtrx,th1,th2);



end