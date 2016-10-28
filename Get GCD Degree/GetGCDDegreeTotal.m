function [t] = GetGCDDegreeTotal(fxy_matrix,gxy_matrix,m,n)
% Calculate the degree of the GCD of two bivariate Power Basis polynomials.
%
% %                 Inputs
%
% fxy_matrix    :
%
% gxy_matrix    :
%
% m
%
% n

% % 
% Calculate the Degree of the GCD

% Initialise some empty vectors for storing during the loop
v_MinimumSingularValue = zeros(min(m,n),1);

v_maxDiagR1 = zeros(min(m,n),1);
v_minDiagR1 = zeros(min(m,n),1);

v_maxRowNormR1 = zeros(min(m,n),1);
v_minRowNormR1 = zeros(min(m,n),1);

% initialise some useful vectors
Data_RowNorm    = [];
Data_DiagNorm   = [];

lower_lim = 1;
upper_lim = min(m,n);


%%
% pad the coefficients of fxy and gxy
% this is equivalent to degree elevating so that f is of degree (m,m), and
% g is of degree (n,n)
fxy_matrix_padd = zeros(m+1,m+1);
gxy_matrix_padd = zeros(n+1,n+1);


[r,c] = size(fxy_matrix);
fxy_matrix_padd(1:r,1:c) = fxy_matrix;

[r,c] = size(gxy_matrix);
gxy_matrix_padd(1:r,1:c) = gxy_matrix;

%%
data = [];
% let k represent the total degree of the common divisor
for k = 1:1:min(m,n)
    
    % Build the partitions of the Sylvester matrix
    T1 = BuildT1_TotalDegree(fxy_matrix_padd,m,n-k);
    T2 = BuildT1_TotalDegree(gxy_matrix_padd,n,m-k);
    
    % Build the sylvester matrix
    Sk = [T1 T2];
    
    % Using QR Decomposition of the sylvester matrix
    [~,R] = qr(Sk);
    
    % Take absolute values.
    R = abs(R);
    
    % Store all of the diagonals of R
    vec_diags = diag(R);
    ks = k.* ones(length(vec_diags),1);
    
    data = [data ; ks vec_diags];
    
    % Get number of rows in R1
    [R1_rows,~] = size(diag(R));
    
    % Obtain R1 the top square of the R matrix.
    R1 = R(1:R1_rows,1:R1_rows);
    
    % Get Norms of each row in the matrix R1
    R1_RowNorm = sqrt(sum(R1.^2,2))./norm(R1);
    
    % Get ONLY the diagonal elements and normalise them.
    R1_DiagNorm = diag(R1)./norm(diag(R1));
    
    % Scatter Plot Data
    ks = k.*ones(size(R1_RowNorm));
    ns = 1:1:size(R1_RowNorm,1);
    
    % Form a triple of k, QR_RowNorm, and index of the value of
    % the row of R1 corresponding to QR_RowNorm].
    % EG.
    %  [1   0.015  1
    %   1   0.156  2
    %   2 ...]
    X = [ks R1_RowNorm ns'];
    Data_RowNorm = [Data_RowNorm; X];
    
    X2 = [ks R1_DiagNorm ns'];
    Data_DiagNorm = [Data_DiagNorm;X2];
    
    
    % Get the maximum diagonal of R1
    v_maxDiagR1(k) = max(diag(R1));
    
    % Get the minimum diagonal of R1
    v_minDiagR1(k) = min(diag(R1));
    
    % Get the maximum row norm
    v_maxRowNormR1(k) = max(R1_RowNorm);
    
    % Get the minimal row norm
    v_minRowNormR1(k) = min(R1_RowNorm);
    
    % Get minimum singular value of S_{k}
    v_MinimumSingularValue(k) = min(svd(Sk));
    
    
end

vRatio_MaxMin_Diags_R1 = v_maxDiagR1./v_minDiagR1;

vRatio_MaxMin_RowNorm_R1 = v_maxRowNormR1./v_minRowNormR1;


if k == 1 
   fprintf('Only one subresultant exists')
end


% if only one subresultant exists.
if lower_lim == upper_lim
    t = GetGCDDegree_OneSubresultant(Sk);
    return;
else
    t = GetGCDDegree_MultipleSubresultants(v_MinimumSingularValue,[1,min(m,n)]);
end
   
% Plot graphs
PlotGraphs()

end
