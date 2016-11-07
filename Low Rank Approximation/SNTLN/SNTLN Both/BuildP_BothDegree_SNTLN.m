function Pk = BuildP_BothDegree_SNTLN(m,m1,m2,n,n1,n2,k,k1,k2,alpha,th1,th2,idx_col)
% Calculate the matrix P where P is the matrix such that a column of the 
% Sylvester subresultant matrix S_{k,k1,k2} can be written as a product of
% P and the column vector of coefficients of f and g.
%
% Inputs
%
% m : Total degree of polynomial f(x,y)
%
% m1 : Degree of f(x,y) with respect to x
%
% m2 : Degree of f(x,y) with respect to y
%
% n : Total degree of polynomial g(x,y)
%
% n1 : Degree of g(x,y) with respect to x
%
% n2 : Degree of g(x,y) with respect to y
%
% k : Total degree of d(x,y)
%
% k1 : Degree of d(x,y) with respect to x
%
% k2 : Degree of d(x,y) with respect to y
%
% alpha : Optimal value of \alpha
%
% th1 : Optimal value of \theta_{1}
%
% th2 : Optimal value of \theta_{2}
%
% idx_col : Index of column removed from S_{k_{1},k_{2}}
%
%
% % Outputs
%
% P : Matrix P

% Get the number of coefficients in polynomial f
% % nCoeff_f = (m1+1).*(m2+1);
nCoeff_fxy = GetNumNonZeros(m1,m2,m);

% Get the number of coefficients in polynomial g(x,y)
% % nCoeff_g = (n1+1).*(n2+1);
nCoeff_gxy = GetNumNonZeros(n1,n2,n);

% Get number of Rows in Sylvester matrix S_{k,k1,k2}
nRows_Skk1k2 = GetNumNonZeros(m1+n1-k1,m2+n2-k2,m+n-k);


% Get number of columns in first partition of the Sylvester subresultant
% matrix.
nCols_T1 = GetNumNonZeros(n1-k1,n2-k2,n-k);

if idx_col <= nCols_T1
    % Optimal column in first partition
    
    % % Build the matrix P
    
    % Build the matrix P1
    P1 = BuildP1_BothDegree_SNTLN(m,m1,m2,n,n1,n2,k,k1,k2,th1,th2,idx_col);
    
    % Build the matrix P2
    P2 = zeros(nRows_Skk1k2,nCoeff_gxy);
    
    % Build the matrix P
    Pk = [P1 P2];
    
else
    % Optimal column in second partition
    
    % Build the matrix P1
    P1 = zeros(nRows_Skk1k2,nCoeff_fxy);
    
    % Build the matrix P2
    % Get the position of the optimal column with respect to T(g)
    opt_col_rel = idx_col - nCols_T1;
    
    P2 = BuildP1_BothDegree_SNTLN(n,n1,n2,m,m1,m2,k,k1,k2,th1,th2,opt_col_rel);
    
    % Build the matrix P.
    Pk = [P1 alpha.*P2];
    
end

end
