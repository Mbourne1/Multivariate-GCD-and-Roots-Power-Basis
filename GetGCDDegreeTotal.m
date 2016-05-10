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

global SETTINGS

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
    
    
    %
    v_maxDiagR1(k) = max(diag(R1));
    
    %
    v_minDiagR1(k) = min(diag(R1));
    
    %
    v_maxRowNormR1(k) = max(R1_RowNorm);
    
    %
    v_minRowNormR1(k) = min(R1_RowNorm);
    
    % Get minimum singular value of S_{k}
    v_MinimumSingularValue(k) = min(svd(Sk));
    
    
end

vRatio_MaxMin_Diags_R1 = v_maxDiagR1./v_minDiagR1;

vRatio_MaxMin_RowNorm_R1 = v_maxRowNormR1./v_minRowNormR1;


if k == 1 
   fprintf('Only one subresultant exists')
end

% first check if only one subresultant exists
[r,~] = size(v_MinimumSingularValue);
if (r == 1)
    fprintf('Only one subresultant exists. Check if near rank deficient \n')
    
    vSingularValues = svd(Sk);
    
    % Plot all singular values of S_{1}(f,g)
    figure('name','GetDegreeTotal - SVD')
    hold on
    plot(log10(vSingularValues),'-s')
    hold off
    
    vDelta_MinSingularValues = abs(diff(log10(vSingularValues)));
    
    max_change = max(abs(diff(log10(vSingularValues))));
    
    if max_change < THRESHOLD
        %
        fprintf('Change in Singular values is not significant \n')
        t = 0;
        
        
    else
        fprintf('Change in Singular values is signficant \n')
        t = 1;
    end
    
    fprintf('Degree of GCD : %i \n',t)
    return
end

switch SETTINGS.PLOT_GRAPHS
    case 'y'
        figure_name = sprintf('%s - R Diagonals',mfilename);
        figure('name',figure_name)
        hold on
        title('Diagonal entries of R matrix where S_{k,k} = QR')
        xlabel('k')
        ylabel('log_{10} diagonal r_{i,i}')
        scatter(data(:,1),log10(data(:,2)))
        hold off
        
        figure_name = sprintf('%s - Minimum Singular Values',mfilename);
        figure('name',figure_name)
        titleString = sprintf('Singular Values');
        title(titleString)
        xlabel('k: total degree')
        ylabel('log_{10} \sigma_{i}')
        hold on
        plot(log10(v_MinimumSingularValue),'-s');
        hold off
        
        % plot all the largest ratios for k = 1,...,min(m,n)
        figure_name = sprintf('%s - QR max:min diagonals',mfilename);
        figure('name',figure_name)
        hold on
        title('Plotting max:min diagonal entries of QR decomposition of S_{k}')
        plot(log10(vRatio_MaxMin_Diags_R1),'-s');
        xlabel('k')
        ylabel('log_{10}')
        hold off
        
        figure_name = sprintf('%s - QR max:min Row Sum',mfilename);
        figure('name',figure_name)
        hold on
        title('Plotting max:min rowsums of QR decomposition of S_{k}')
        plot(log10(vRatio_MaxMin_RowNorm_R1),'-s');
        xlabel('k')
        ylabel('log_{10}')
        hold off
    case 'n'
    otherwise
        error('plot_graphs is either y or n')
end

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
