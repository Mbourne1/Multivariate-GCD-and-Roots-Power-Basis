function Pt = BuildP_TotalDegree_STLN(m,n,idx_col,k)
% BuildPt(m,m1,m2,n,n1,n2,opt_col,t1,t2)
%
% Build the matrix P_{t}, such that the matrix vector product P*[f;g] gives
% the column c_{t}.
%
% P_{t} * [f;g] = c_{t}
%
% Inputs
%
% m : Degree of polynomial f(x,y) 
%
% n : Degree of polynomial g(x,y)
%
% idx_col :
%
% k : Degree of polynomial d(x,y) 

% Get the number of non-zero coefficients in polynomial f
nCoeff_f = nchoosek(m+2,2);

% Get the number of non-zero coefficients in polynomial g
nCoeff_g = nchoosek(n+2,2);

% Number of columns in T1 of the sylvester matrix
nColumnsT1 = nchoosek(n-k+2,2);
nColumnsT2 = nchoosek(m-k+2,2);

nRowsT1 = nchoosek(m+n-k+2,2);

if idx_col <= nColumnsT1
    % Optimal column in first partition
    
    % % Build the matrix P
    
    % Build the matrix P1
    P1 = BuildP1_TotalDegree_STLN(m,n,idx_col,k);
    
    % Build the matrix P2
    P2 = zeros(nRowsT1,nCoeff_g);
    
    
    
else
    % Optimal column in second partition
    
    % Build the matrix P1
    P1 = zeros(nRowsT1,nCoeff_f);
    
    % Build the matrix P2
    % Get the position of the optimal column with respect to T(g)
    idx_col_rel = idx_col - nColumnsT1;
    
    P2 = BuildP1_TotalDegree_STLN(n,m,idx_col_rel,k);
    
    
end

Pt = [P1 P2];

end