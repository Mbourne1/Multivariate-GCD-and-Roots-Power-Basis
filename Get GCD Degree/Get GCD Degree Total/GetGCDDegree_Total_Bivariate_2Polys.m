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


% let k represent the total degree of the common divisor d_{k}(x,y)
for i = 1 : 1 : nSubresultants
    
    k = lowerLimit + (i-1) ;
    

    % Build the partitions of the Sylvester matrix
    T_f = BuildT1_Total_Bivariate(fxy_matrix_padd, m, n-k);
    T_g = BuildT1_Total_Bivariate(gxy_matrix_padd, n, m-k);


    % Build the sylvester matrix
    Sk = [T_f T_g];
        
        
 
    
    % Using QR Decomposition of the sylvester matrix
    [~,R] = qr(Sk);
    
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
    arr_SingularValues{i} = svd(Sk);
    
end



% R1 Row Norms
% R1 Row Diagonals
% Singular Values
% Residuals
global SETTINGS

switch SETTINGS.RANK_REVEALING_METRIC
    
    case 'Singular Values'
        
        % Initialise a vector to store minimum singular values
        vMinimumSingularValue = zeros(nSubresultants,1);
        
        for i = 1:1:nSubresultants
            vMinimumSingularValue(i) = min(arr_SingularValues{i});
        end
        
        if (SETTINGS.PLOT_GRAPHS)
            plotMinimumSingularValues_degreeTotal(vMinimumSingularValue, limits);
            plotSingularValues_degreeTotal(arr_SingularValues, limits);
        end
        
        metric = vMinimumSingularValue;
        
    case 'R1 Row Norms'
        % Get max/min row norms
        
        % Initialise vectors to store maximum row norm and minimum row norm
        vMaxRowNormR1 = zeros(nSubresultants,1);
        vMinRowNormR1 = zeros(nSubresultants,1);
        
        for i = 1:1:nSubresultants
            % Get the maximum row norm
            vMaxRowNormR1(i) = max(arr_R1_RowNorm{i});
            
            % Get the minimal row norm
            vMinRowNormR1(i) = min(arr_R1_RowNorm{i});
        end
        
        vRatio_MaxMin_RowNorm_R1 = vMinRowNormR1./vMaxRowNormR1;
        
        metric = vRatio_MaxMin_RowNorm_R1;
        
        if (SETTINGS.PLOT_GRAPHS)
            plotRowNorm_degreeTotal(arr_R1_RowNorm, limits)
            plotMaxMinRowNorm_degreeTotal(vRatio_MaxMin_RowNorm_R1, limits);
        end
        
    case 'R1 Row Diagonals'
        
        % Initialise vectors
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
    t = GetGCDDegree_OneSubresultant(Sk);
    return;
else
    t = GetGCDDegree_MultipleSubresultants(metric, limits );
end

% Plot graphs
%PlotGraphs()

end
