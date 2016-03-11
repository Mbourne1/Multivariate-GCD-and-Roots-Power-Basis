function [t1,t2] = getGCDDegree(fxy_matrix,gxy_matrix,m,n)
% Calculate the degree of the GCD of two bivariate Power Basis polynomials.

% %                 Inputs

% fxy_matrix

% gxy_matrix

% Opt_theta_1

% Opt_theta_2

%%
% Get degrees of fxy 
[rows,cols] = size(fxy_matrix);

% Get degree of f in terms of x
m1 = rows - 1;

% Get degree of f in terms of y
m2 = cols - 1;

% Get degree of gxy
[rows,cols] = size(gxy_matrix);

% Get degree of g in terms of x
n1 = rows - 1;

% Get degree of g in terms of y
n2 = cols - 1;

Y=[];

m = m1 + m2
n = n1 + n2


%% Calculate the Degree of the GCD

% Initialise some empty vectors for storing during the loop
sing_val_vec = zeros(min(m,n),1);
t_val_vec1 = zeros(min(m,n),2);
t_val_vec2 = zeros(min(m,n),2);
t_val_vec3 = zeros(min(m,n),2);

% initialise some useful vectors
Data_RowNorm    = [];
Data_DiagNorm   = [];

    
% let k represent the 
for k = 1:1:min(m,n)
    min_sing_val_vec = [];
    temp_vec = [];
    count = 1;
    ratio_maxmin_diag_vec = [];
    ratio_maxmin_rowsum_vec = [];
    
    % fix the degree with respect to x, then cycle through all degrees for
    % y. 
    for k1 = 0:1:k
        k2 = k-k1;
        if (k1 > m1 || k1 > n1) || (k2 > m2 || k2 > n2)
        else
     

            T1 = BuildT1(fxy_matrix,n1,n2,k1,k2);
            T2 = BuildT1(gxy_matrix,m1,m2,k1,k2);
            Sk = [T1 T2];

            
            fprintf('max and min entry for Subresultant S%i,%i \n ', k1,k2)
            fprintf('Max : %0.5e \n', max(max(abs(Sk))) );
            fprintf('Max : %0.5e \n', min(min(abs(Sk))) );
                    
            
            [rows,cols] = size(Sk);
            
            % Edit 24/07/2015
            % Using QR Decomposition of the sylvester matrix
            [~,R] = qr(Sk);
            
            % Take absolute values.
            R = abs(R);
            
            % Get number of rows in R1
            [R1_rows,~] = size(diag(R));
            
            % Obtain R1 the top square of the R matrix.
            R1 = R(1:R1_rows,1:R1_rows);
            
            % Get Norms of each row in the matrix R1
            R1_RowNorm = sqrt(sum(R1.^2,2))./norm(R1);
            
            % Get ONLY the diagonal elements and normalise them.
            R1_DiagNorm = diag(R1)./norm(diag(R1));
            
            % Scatter Plot Data
            ks = count.*ones(size(R1_RowNorm));
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
            

            
            % add the minimum singular value S_{k1,k2} to the vector of
            % singular values for all k = k1+k2
            min_sing_val_vec(count) = min(svd(Sk));
            temp_vec(count,:) = [k1 k2];
            count = count + 1;
            
        end
    end
    figure(999)
    outputString = sprintf('Minimal singular value for all subresultants of degree k = %i', k);
    title(outputString)
    hold on
    plot(log10(min_sing_val_vec));
    hold off
%     
%     figure(998)
%     hold on
%     scatter(Data_RowNorm(:,1), Data_RowNorm(:,2))
    
    
    
    try
    % of all the diagonal entries in R1 for the given k1+k2 = k (where k is the total degree), get the
    % largest ratio
    max(ratio_maxmin_diag_vec)
    
    [ratio_maxmin_diag(k),index1] = max(ratio_maxmin_diag_vec)
    t_val_vec1(k,:) = temp_vec(index1,:);
    
    [ratio_maxmin_rowsum(k),index2] = max(ratio_maxmin_rowsum_vec);
    t_val_vec2(k,:) = temp_vec(index2,:);
    
    % of all the minimal singular values for the given k1+k2 = k, stored in
    % the vector min_sing_val_vec, get the minimal singular value and its
    % index.
    [sing_val_vec(k),index] = min(min_sing_val_vec);
    
    % get the values of k1 and k2 which gave the minimal value
    t_val_vec3(k,:) = temp_vec(index,:);
    catch
    fprintf('No Values')    
    end
end
        
% plot all the largest ratios for k = 1,...,min(m,n)
figure(1)
hold on
title('Plotting max:min diagonal entries of QR decomposition of S_{k}')
plot(log10(ratio_maxmin_diag),'-s');
xlabel('k')
ylabel('log_{10}')
hold off

figure(2)
hold on
title('Plotting max:min rowsums of QR decomposition of S_{k}')
plot(log10(ratio_maxmin_rowsum),'-s');
xlabel('k')
ylabel('log_{10}')
hold off

% plot all the minimum singular values for k = 1,...,min(m,n)
figure(3)
hold on
title('Singular values of each Subresultant S_{k}')
plot(log10(sing_val_vec),'-s')
xlabel('k')
ylabel('Smallest Singular Value')
hold off




[~,maxindex] = max(diff(log10(sing_val_vec)));
degree_calc = maxindex;

[tval] = t_val_vec3(degree_calc,:);
t1 = tval(1);
t2 = tval(2);

    
        
end
