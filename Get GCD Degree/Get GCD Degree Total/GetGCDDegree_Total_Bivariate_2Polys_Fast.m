function [t] = GetGCDDegree_Total_Bivariate_2Polys(fxy, gxy, m, n, limits)
% Calculate the degree of the GCD of two bivariate power basis polynomials.
%
% % Inputs
%
% [fxy, gxy] : Coefficients of polynomial f(x,y) and g(x,y)
%
% [m, n] : Total degree of polynomial f(x,y) and g(x,y)
%
% % Outputs
%
% t : Total degree of the GCD d(x,y)

% %
% Calculate the Degree of the GCD

lowerLimit = limits(1);
upperLimit = limits(2);


nSubresultants = upperLimit - lowerLimit + 1;

% Initialise some cell arrays
arr_R1_RowNorm = cell(nSubresultants,1);
arr_R1_Diag = cell(nSubresultants,1);
arr_SingularValues = cell(nSubresultants,1);




%%
% pad the coefficients of fxy and gxy
% this is equivalent to degree elevating so that f is of degree (m,m), and
% g is of degree (n,n)
fxy_matrix_padd = zeros(m+1,m+1);
gxy_matrix_padd = zeros(n+1,n+1);

[r,c] = size(fxy);
fxy_matrix_padd(1:r,1:c) = fxy;

[r,c] = size(gxy);
gxy_matrix_padd(1:r,1:c) = gxy;

arr_Sk = cell(nSubresultants,1);

% let k represent the total degree of the common divisor d_{k}(x,y)
for i = 1 : 1 : nSubresultants
    
    k = lowerLimit + (i-1) ;
    
    if (i == 1)
        
    % Build the partitions of the Sylvester matrix
    T_f = BuildT1_Total_Bivariate(fxy_matrix_padd, m, n-k);
    T_g = BuildT1_Total_Bivariate(gxy_matrix_padd, n, m-k);
    
    % Build the sylvester matrix
    arr_Sk{i} = [T_f T_g];
    
    else
        
        nCols_T1_prev = nchoosek(n-(k-1)+2,2);
        nCols_T1_Removed = n-(k-1) + 1;
        idx_firstColumnRemoved_T1 = nCols_T1_prev - nCols_T1_Removed + 1;
        idx_lastColumnRemoved_T1 = nCols_T1_prev;
        
        % Get vector of column index to be removed from T(f)
        vIndexColsRemoved_T1 = [idx_firstColumnRemoved_T1 : 1 : idx_lastColumnRemoved_T1];
        
        nCols_T2_prev = nchoosek(m-(k-1)+2,2); 
        nCols_Sk_prev = nCols_T1_prev + nCols_T2_prev;
       
        nCols_T2_Removed = m - (k-1) + 1;
        idx_firstColumnRemoved_T2 = nCols_Sk_prev - nCols_T2_Removed + 1;
        idx_lastColumnRemoved_T2 = nCols_Sk_prev;
        
        % Get vector of column index to be removed from T(g)
        vIndexColsRemoved_T2 = [idx_firstColumnRemoved_T2 : 1 : idx_lastColumnRemoved_T2];
        
        nRows_prev = nchoosek(m+n-(k-1)+2,2);
        nRows_Removed = m+n-(k-1)+1;
        idx_firstRowRemoved = nRows_prev - nRows_Removed +1;
        idx_lastRowRemoved = nRows_prev;
        
        
        Sk_prev = arr_Sk{i-1};
        
        % Remove n-k+1 columns from the end of T_{n-k}(f(x,y)) to get
        % T_{n-k-1}(f(x,y))
        % Remove m-k+1 columns from the end of T_{m-k}(g(x,y)) to get
        % T_{m-k+1}(f(x,y))
        Sk_prev(:,[vIndexColsRemoved_T1 vIndexColsRemoved_T2]) = [];
        
        % Remove m+n-k+1 rows from the bottom of both partitions
        Sk_prev(idx_firstRowRemoved:idx_lastRowRemoved,:) = [];
        
        arr_Sk{i} = Sk_prev;
              
        
    end
    
    % Using QR Decomposition of the sylvester matrix
    [R] = qr(arr_Sk{i},0);
    
    % Take absolute values.
    R = abs(R);
    
    % Get number of rows in R1
    [R1_rows,~] = size(diag(R));
    
    % Obtain R1 the top square of the R matrix.
    R1 = R(1:R1_rows,1:R1_rows);
    
    % Get Norms of each row in the matrix R1
    arr_R1_RowNorm{i} = sqrt(sum(R1.^2,2))./norm(R1);
    
    % Get ONLY the diagonal elements and normalise them.
    arr_R1_Diag{i} = diag(R1);
    
    % Get singular values of S_{k}
    arr_SingularValues{i} = svd(arr_Sk{i});
    
end



% R1 Row Norms
% R1 Row Diagonals
% Singular Values
% Residuals
global SETTINGS

switch SETTINGS.RANK_REVEALING_METRIC
    
    case 'Singular Values'
        
        vMinimumSingularValue = zeros(nSubresultants,1);
        
        for i = 1:1:nSubresultants
            vMinimumSingularValue(i) = min(arr_SingularValues{i});
        end
        
        if(SETTINGS.PLOT_GRAPHS)
            plotMinimumSingularValues_degreeTotal(vMinimumSingularValue, limits);
            plotSingularValues_degreeTotal(arr_SingularValues, limits);
        end
        metric = vMinimumSingularValue;
        
    case 'R1 Row Norms'
        % Get max/min row norms
        
        vMaxRowNormR1 = zeros(nSubresultants,1);
        vMinRowNormR1 = zeros(nSubresultants,1);

        for i = 1:1:nSubresultants
            % Get the maximum row norm
            vMaxRowNormR1(i) = max(arr_R1_RowNorm{i});
            
            % Get the minimal row norm
            vMinRowNormR1(i) = min(arr_R1_RowNorm{i});
        end
        
        vRatio_MaxMin_RowNorm_R1 = vMaxRowNormR1./vMinRowNormR1;
        metric = vRatio_MaxMin_RowNorm_R1;
        
        if(SETTINGS.PLOT_GRAPHS)
            plotRowNorm_degreeTotal(arr_R1_RowNorm, limits)
            plotMaxMinRowNorm_degreeTotal(vRatio_MaxMin_RowNorm_R1, limits);
        end
        
    case 'R1 Row Diagonals'
        
        v_maxDiagR1 = zeros(nSubresultants,1);
        v_minDiagR1 = zeros(nSubresultants,1);
        
        for i = 1:1:nSubresultants
            % Get the maximum diagonal of R1
            v_maxDiagR1(i) = max(abs(arr_R1_Diag{i}));
            
            % Get the minimum diagonal of R1
            v_minDiagR1(i) = min(abs(arr_R1_Diag{i}));
        end
        
        
        % Get max/min diagonal entries
        vRatio_MaxMin_Diags_R1 = v_minDiagR1./v_maxDiagR1;
        metric = vRatio_MaxMin_Diags_R1;
        
        if(SETTINGS.PLOT_GRAPHS)
            plotRowDiag_degreeTotal(arr_R1_Diag, limits);
            plotMaxMinRowDiag_degreeTotal(vRatio_MaxMin_Diags_R1, limits);
        end
        
    case 'Residuals'
        error('Not Developed')
end


if i == 1
    fprintf('Only one subresultant exists')
end


% if only one subresultant exists.
if lowerLimit == upperLimit
    t = GetGCDDegree_OneSubresultant(arr_Sk);
    return;
else
    t = GetGCDDegree_MultipleSubresultants(metric, limits );
end

% Plot graphs
%PlotGraphs()

end
