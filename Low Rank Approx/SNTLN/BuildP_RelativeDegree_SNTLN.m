function P = BuildP_RelativeDegree_SNTLN(m1,m2,n1,n2,alpha,th1,th2,idx_opt_col,k1,k2,nCols_T1)
% Calculate the matrix DP where P is the matrix such that c = P[f;g]
%
% Inputs
%
% m1 : Degree of f(x,y) with respect to x
%
% m2 : Degree of f(x,y) with respect to y
%
% n1 : Degree of g(x,y) with respect to x
%
% n2 : Degree of g(x,y) with respect to y
%
% alpha : Optimal value of \alpha
%
% th1 : Optimal value of \theta_{1}
%
% th2 : Optimal value of \theta_{2}
%
% idx_opt_col : Index of column removed from S_{k_{1},k_{2}}
%
% k1 : Degree of d(x,y) with respect to x
%
% k2 : Degree of d(x,y) with respect to y
%
% num_cols_T1 : Number of

% Get the number of coefficients in polynomial f
nCoeff_f = (m1+1).*(m2+1);

% Get the number of coefficients in polynomial g
nCoeff_g = (n1+1).*(n2+1);

nRows = (m1+n1-k1+1)*(m2+n2-k2+1);

if idx_opt_col <= nCols_T1
    % Optimal column in first partition
    
    % % Build the matrix P
    
    % Build the matrix P1
    P1 = BuildP1_RelativeDegree_SNTLN(m1,m2,n1,n2,th1,th2,idx_opt_col,k1,k2);
    
    % Build the matrix P2
    P2 = zeros(nRows,nCoeff_g);
    
    % Build the matrix P
    P = [P1 P2];
    
else
    % Optimal column in second partition
    
    % Build the matrix P1
    P1 = zeros(nRows,nCoeff_f);
    
    % Build the matrix P2
    % Get the position of the optimal column with respect to T(g)
    opt_col_rel = idx_opt_col - nCols_T1;
    P2 = BuildP1_RelativeDegree_SNTLN(n1,n2,m1,m2,th1,th2,opt_col_rel,k1,k2);
    
    % Build the matrix P.
    P = [P1 alpha.*P2];
    
end

end
