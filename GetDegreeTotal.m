function [t] = GetDegreeTotal(fxy_matrix,gxy_matrix,m,n,...
    lambda,mu,...
    opt_alpha, opt_theta_1,opt_theta_2)
% Calculate the degree of the GCD of two bivariate Power Basis polynomials.
%
% %                 Inputs
%
% fxy_matrix    :
%
% gxy_matrix    :
%
% opt_alpha     :
%
% Opt_theta_1   :
%
% Opt_theta_2   :
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global PLOT_GRAPHS
global THRESHOLD % USED ON LINE


% Get degrees of f(x,y)
[rows,cols] = size(fxy_matrix);

% Get degree of f(x,y) in terms of x
m1 = rows - 1;

% Get degree of f(x,y) in terms of y
m2 = cols - 1;

% Get degree of g(x,y)
[rows,cols] = size(gxy_matrix);

% Get degree of g(x,y) in terms of x
n1 = rows - 1;

% Get degree of g(x,y) in terms of y
n2 = cols - 1;

%% 
% Preprocess the polynomials

% Normalise Polynomial f(x,y) by geometric mean
fxy_matrix_n = fxy_matrix ./ lambda;

% Normalise polynomial g(x,y) by geometric mean
gxy_matrix_n = gxy_matrix ./ mu;

% Preprocess Polynomial f(x,y) to obtain f(w,w)
th1_mat = diag(opt_theta_1.^(0:1:m1));
th2_mat = diag(opt_theta_2.^(0:1:m2));
fww_matrix = th1_mat * fxy_matrix_n * th2_mat;

% Preprocess polynomial g(x,y) to obtain g(w,w)
th1_mat = diag(opt_theta_1.^(0:1:n1));
th2_mat = diag(opt_theta_2.^(0:1:n2));
gww_matrix = th1_mat * gxy_matrix_n * th2_mat;


%% 
% Calculate the Degree of the GCD

% Initialise some empty vectors for storing during the loop
sing_val_vec = zeros(min(m,n),1);

% initialise some useful vectors
Data_RowNorm    = [];
Data_DiagNorm   = [];

ratio_maxmin_diag_vec = [];
ratio_maxmin_rowsum_vec = [];
vMin_sing_val_vec = [];

%%
% pad the coefficients of fxy and gxy
% this is equivalent to degree elevating so that f is of degree (m,m), and
% g is of degree (n,n)
fww_matrix_padd = zeros(m+1,m+1);
gww_matrix_padd = zeros(n+1,n+1);

[r,c] = size(fww_matrix);
fww_matrix_padd(1:r,1:c) = fww_matrix;

[r,c] = size(gww_matrix);
gww_matrix_padd(1:r,1:c) = gww_matrix;

%%
data = [];
% let k represent the total degree of the common divisor
for k = 1:1:min(m,n)
    
    % Build the partitions of the Sylvester matrix
    T1 = BuildT1_totaldegree(fww_matrix_padd,m,n,k);
    T2 = BuildT1_totaldegree(gww_matrix_padd,n,m,k);
    
    % Build the sylvester matrix
    Sk = [T1 opt_alpha.*T2];
    
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
    
    % Get ratio of max diag elem of S_{k} to min diag elemente of S_{k}
    ratio_maxmin_diag_vec = [ratio_maxmin_diag_vec ; max(diag(R1))./min(diag(R1))];
    ratio_maxmin_rowsum_vec = [ratio_maxmin_rowsum_vec ; max(R1_RowNorm)./min(R1_RowNorm)];
    
    % Get the minimum singular value
    vMin_sing_val_vec = [vMin_sing_val_vec; min(svd(Sk))];
    
    
    
    
end





% first check if only one subresultant exists
[r,~] = size(vMin_sing_val_vec);
if (r == 1)
    fprintf('Only one subresultant exists. Check if near rank deficient \n')
    
    vSingularValues = svd(Sk)
    
    % Plot all singular values of S_{1}(f,g)
    figure('name','GetDegreeTotal - SVD')
    hold on
    plot(log10(vSingularValues),'-s')
    hold off
    
    vDelta_MinSingularValues = abs(diff(log10(vSingularValues)))
    
    max_change = max(abs(diff(log10(vSingularValues))))
    
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

switch PLOT_GRAPHS
    case 'y'
        figure('name','Get Degree Total : R Diagonals')
        hold on
        title('Diagonal entries of R matrix where S_{k,k} = QR')
        xlabel('k')
        ylabel('log_{10} diagonal r_{i,i}')
        scatter(data(:,1),log10(data(:,2)))
        hold off
        
        figure('name','Minimum Singular values')
        titleString = sprintf(['Minimal singular values of each subresultant S_{k} \n'...
            'alpha  = %0.5e ' ...
            'theta1 = %0.5e ' ...
            'theta2 = %0.5e '], ...
            opt_alpha , opt_theta_1 , opt_theta_2 );
        title(titleString)
        xlabel('k: total degree')
        ylabel('log_{10} \sigma_{i}')
        hold on
        plot(log10(vMin_sing_val_vec),'-s');
        hold off
        
        % plot all the largest ratios for k = 1,...,min(m,n)
        figure('name','Get Degree Total : R max:min diagonals')
        hold on
        title('Plotting max:min diagonal entries of QR decomposition of S_{k}')
        plot(log10(ratio_maxmin_diag_vec),'-s');
        xlabel('k')
        ylabel('log_{10}')
        hold off
        
        figure('name','GetDegreeTotal : QR max:min row sum')
        hold on
        title('Plotting max:min rowsums of QR decomposition of S_{k}')
        plot(log10(ratio_maxmin_rowsum_vec),'-s');
        xlabel('k')
        ylabel('log_{10}')
        hold off
    case 'n'
    otherwise
        error('plot_graphs is either y or n')
end

% Calculate the total degree

fprintf('-----------------------------------------------------------------\n')
[svd_val,svd_maxindex] = max(diff(log10(vMin_sing_val_vec)));
fprintf('Total Degree Calculated By Minimum Singular Values: %i \n',svd_maxindex);

[rowdiag_val,rowdiag_maxindex] = min(diff(log10(ratio_maxmin_diag_vec)));
fprintf('Total Degree Calculated By Max:Min Row Diags: %i \n',rowdiag_maxindex);

[rowsum_val,rowsum_maxindex] = min(diff(log10(ratio_maxmin_rowsum_vec)));
fprintf('Total Degree Calculated By Max:Min Row Sums: %i \n',rowsum_maxindex);

fprintf('-----------------------------------------------------------------\n')


val = rowdiag_val;
index = rowdiag_maxindex;

% check if the maximum change is significant

if abs(val) < THRESHOLD
    
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
