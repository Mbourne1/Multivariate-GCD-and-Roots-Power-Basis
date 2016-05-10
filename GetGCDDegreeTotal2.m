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

global SETTINGS % USED ON LINE



lower_lim = limits_t(1);
upper_lim = limits_t(2);

if (lower_lim == upper_lim)
    t = lower_lim;
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

n_Subresultants = upper_lim - lower_lim + 1;

v_maxDiagR1 = zeros(n_Subresultants,1);
v_minDiagR1 = zeros(n_Subresultants,1);
v_maxRowNormR1 = zeros(n_Subresultants,1);
v_minRowNormR1 = zeros(n_Subresultants,1);
v_MinimumSingularValue = zeros(n_Subresultants,1);

% %
data = [];


% let k represent the total degree of the common divisor
for k = lower_lim:1:upper_lim
    
    % Get i, an index for any storage vectors which will always start at 1.
    i = k - lower_lim + 1;
    
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

% Get ratio of max to min entry on the diagonal of R1_{i} for each
% subresultant S_{i} = QR_{i}.
vRatio_MaxMin_Diags_R1 = v_maxDiagR1./v_minDiagR1;

% Get the ratio of max to min row norm in R1_{i} for each subresultant S_{i} = QR_{i}
vRatio_MaxMin_RowNorm_R1 = v_maxRowNormR1./v_minRowNormR1;

% if only one subresultant exists.
if lower_lim == upper_lim
    GetRank_One_Subresultant(vSingularValues);
else
end
   

plotgraphs()


% Calculate the total degree

fprintf('-----------------------------------------------------------------\n')
[svd_val,svd_maxindex] = max(diff(log10(v_MinimumSingularValue)));
fprintf('Total Degree Calculated By Minimum Singular Values: %i \n',svd_maxindex);

[rowdiag_val,rowdiag_maxindex] = min(diff(log10(vRatio_MaxMin_Diags_R1)));
fprintf('Total Degree Calculated By Max:Min Row Diags: %i \n',rowdiag_maxindex);

[rowsum_val,rowsum_maxindex] = min(diff(log10(vRatio_MaxMin_RowNorm_R1)));
fprintf('Total Degree Calculated By Max:Min Row Sums: %i \n',rowsum_maxindex);

fprintf('-----------------------------------------------------------------\n')


val = rowdiag_val;
index = rowdiag_maxindex;

% check if the maximum change is significant

if abs(val) < SETTINGS.THRESHOLD
    
    % Not significant
    t = min(m,n);
    fprintf('No significant change in minimal row diagaonl between subresultants \n')
    fprintf('All Subresultants are either full rank or rank deficient. \n')
    fprintf('Degree of GCD either zero or min(m,n)\n')
    
    
else
    % Change is significant
    t = index;
end


end



function t = GetRank_One_Subresultant(vMinSingVal)
% Given the vector
% Get the rank, where only one subresultant exists.

global SETTINGS

% Only one subresultant
fprintf('Only one subresultant exists. \n')
switch SETTINGS.PLOT_GRAPHS
    case 'y'
        figure('name','GetDegree - One Subresultant - Singular Values of S1')
        hold on
        title('Singular values of S_{1}')
        plot(log10(vMinSingVal))
        hold off
    case 'n'
end

% Get the differences between singular values of S_{1}
vDiffSingularValues = diff(log10(vMinSingVal));

% Get the index of the largest change in singular values of S_{1}
[deltaSingularValues, ~] = max(vDiffSingularValues);

% If the change is smaller than the predefined threshold value, then plot
% is considered 'flat'.
if deltaSingularValues < SETTINGS.THRESHOLD
    
    % The subresultant is of full rank, in which case t = 0
    t = 0;
    fprintf('The only Subresultant S_{1} appears to be of NonSingular. \n');
    return
    
else % val > threshold
    
    % The subresultant S_{1} is rank deficient, in which case t = 1
    t = 1;
    fprintf('The only Subresultant S_{1} appears to be Singular \n');
    return
    
end
end