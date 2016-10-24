function [t] = GetGCDDegreeTotal2(fxy_matrix,gxy_matrix,m,n,limits_t)
% Calculate the degree of the GCD of two bivariate Power Basis polynomials.
%
% %                 Inputs
%
% fxy_matrix    :
%
% gxy_matrix    :
%
%

% Set upper and lower bound for total degree t.
lower_lim = limits_t(1);
upper_lim = limits_t(2);

% Set upper and lower limit of the degree of the gcd.
lower_lim_comp = 1;
upper_lim_comp = min(m,n);

limits_t_comp = [lower_lim_comp upper_lim_comp];

% if the upper and lower bound are equal, and not equal to one, then set
% the total degree to lower bound. If bound = 1, then it is possible for
% the polynomials to be coprime.
if (lower_lim_comp == upper_lim_comp && lower_lim_comp ~= 1)
    t = lower_lim_comp;
    return;
end



% initialise some useful vectors
Data_RowNorm    = [];
Data_DiagNorm   = [];


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

% Get the number of subresultants
n_Subresultants = upper_lim_comp - lower_lim_comp + 1;

% Initialise some vectors
v_maxDiagR1 = zeros(n_Subresultants,1);
v_minDiagR1 = zeros(n_Subresultants,1);
v_maxRowNormR1 = zeros(n_Subresultants,1);
v_minRowNormR1 = zeros(n_Subresultants,1);
v_MinimumSingularValue = zeros(n_Subresultants,1);

% %
data = [];


% let k represent the total degree of the common divisor
for k = lower_lim_comp:1:upper_lim_comp
    
    % Get i, an index for any storage vectors which will always start at 1.
    i = k - lower_lim_comp + 1;
    
    % Build the partitions T1 and T2 of the Sylvester matrix
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
    
    % Form a triple of [ks, the value of QR_RowNorm, and the index of the value of
    % the row of R1 corresponding to QR_RowNorm].
    % EG.
    %  [1   0.015  1
    %   1   0.156  2
    %   2 ...]
    X = [ks R1_RowNorm ns'];
    Data_RowNorm = [Data_RowNorm; X];
    
    X2 = [ks R1_DiagNorm ns'];
    Data_DiagNorm = [Data_DiagNorm;X2];
    
    % Get maximum entry on the diagonal of R1
    v_maxDiagR1(i) = max(diag(R1));
    
    % Get minimum entry on the diagonal of R1
    v_minDiagR1(i) = min(diag(R1));
    
    % Get maximum Row Norm of R1
    v_maxRowNormR1(i) = max(R1_RowNorm);
    
    % Get minimum Row Norm of R1
    v_minRowNormR1(i) = min(R1_RowNorm);
    
    % Get minimum singular value of S_{k}
    vSingularValues = svd(Sk);
    v_MinimumSingularValue(i) = min(vSingularValues);
   
end


% if only one subresultant exists.
if lower_lim_comp == upper_lim_comp
    t = GetGCDDegree_OneSubresultant(Sk);
    return;
else
    t = GetGCDDegree_MultipleSubresultants(v_MinimumSingularValue,limits_t_comp);
end
   
% Plot graphs
PlotGraphs_DegreeTotal()




end


