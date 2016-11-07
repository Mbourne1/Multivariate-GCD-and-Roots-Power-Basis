function P = BuildP_TotalDegree_SNTLN(m, n, k, alpha, th1, th2, idx_col)
% Calculate the matrix P where P is the matrix such that a column of the 
% Sylvester subresultant matrix S_{k,k1,k2} can be written as a product of
% P and the column vector of coefficients of f and g.
%
% Inputs
%
% m : Total degree of polynomial f(x,y)
%
% n : Total degree of polynomial g(x,y)
%
% k : Total degree of polynomial d(x,y)
%
% alpha : Optimal value of \alpha
%
% th1 : Optimal value of \theta_{1}
%
% th2 : Optimal value of \theta_{2}
%
% idx_col : Index of column removed from S_{k_{1},k_{2}}
% 
% % Outputs
%
% P : Matrix P

% Get the number of coefficients in polynomial f
nCoeff_fxy = nchoosek(m+2,2);

% Get the number of coefficients in polynomial g(x,y)
% % nCoeff_g = (n1+1).*(n2+1);
nCoeff_gxy = nchoosek(n+2,2);

% Get number of Rows in Sylvester matrix S_{k,k1,k2}
nRows_Skk1k2 = nchoosek(m+n-k+2,2);


% Get number of columns in first partition of the Sylvester subresultant
% matrix.
nCols_T1 = nchoosek(n-k+2,2);
nCols_T2 = nchoosek(m-k+2,2);


if idx_col <= nCols_T1
    % Optimal column in first partition
    
    % % Build the matrix P
    
    % Build the matrix P1
    P1 = BuildP1_TotalDegree_SNTLN(m,n,th1,th2,idx_col,k);
    
    % Build the matrix P2
    P2 = zeros(nRows_Skk1k2,nCoeff_gxy);
    
    % Build the matrix P
    P = [P1 P2];
    
else
    % Optimal column in second partition
    
    % Build the matrix P1
    P1 = zeros(nRows_Skk1k2,nCoeff_fxy);
    
    % Build the matrix P2
    % Get the position of the optimal column with respect to T(g)
    opt_col_rel = idx_col - nCols_T1;
    
    P2 = BuildP1_TotalDegree_SNTLN(n,m,th1,th2,opt_col_rel,k);
    
    % Build the matrix P.
    P = [P1 alpha.*P2];
    
end

end