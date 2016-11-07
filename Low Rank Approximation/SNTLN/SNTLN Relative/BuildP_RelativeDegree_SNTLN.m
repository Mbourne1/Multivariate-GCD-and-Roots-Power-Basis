function P = BuildP_RelativeDegree_SNTLN(m1,m2,n1,n2,k1,k2,alpha,th1,th2,idx_col)
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
% % Outputs.
%
% P : Matrix P

% Get the number of coefficients in polynomial f(x,y)
nCoeff_f = (m1+1).*(m2+1);

% Get the number of coefficients in polynomial g(x,y)
nCoeff_g = (n1+1).*(n2+1);

% Get the number of cols in the first partition of the Sylvester
% subresultant matrix S_{k1,k2}
nCols_T1 = (n1-k1+1) * (n2-k2+1);

nRowsSylvester = (m1+n1-k1+1)*(m2+n2-k2+1);

if idx_col <= nCols_T1
    % Optimal column in first partition
    
    % % Build the matrix P
    
    % Build the matrix P1
    P1 = BuildP1_RelativeDegree_SNTLN(m1,m2,n1,n2,k1,k2,th1,th2,idx_col);
    
    % Build the matrix P2
    P2 = zeros(nRowsSylvester,nCoeff_g);
    
    % Build the matrix P
    P = [P1 P2];
    
else
    % Optimal column in second partition
    
    % Build the matrix P1
    P1 = zeros(nRowsSylvester,nCoeff_f);
    
    % Build the matrix P2
    % Get the position of the optimal column with respect to T(g)
    opt_col_rel = idx_col - nCols_T1;
    P2 = BuildP1_RelativeDegree_SNTLN(n1,n2,m1,m2,k1,k2,th1,th2,opt_col_rel);
    
    % Build the matrix P.
    P = [P1 alpha.*P2];
    
end

end
