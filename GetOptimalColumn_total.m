function [opt_col] = GetOptimalColumn_total(fxy_matrix,gxy_matrix,m,n,t)

% % Build the partitions of the Sylvester matrix S_{t}

% Build the first partition containing coefficients of fxy
T1 = BuildT1_totaldegree(fxy_matrix,m,n,t);

% Build the second partition containing coefficients of gxy
T2 = BuildT1_totaldegree(gxy_matrix,n,m,t);

% Concatenate the two partitions
St = [T1 T2];

% From the given subresultant find the optimal column for removal.
[~,cols_T1] = size(St);

% Take the QR decomposition of the Subresultant
[Qk,Rk] = qr(St);


residuals_QR = zeros(cols_T1,1);
for k=1:1:cols_T1
    Sk_temp = St;
    % Rem
    ck = Sk_temp(:,k);
    [Q,~] = qrdelete(Qk,Rk,k);
    cd = Q'*ck;
    d = cd(n+1:end,:);
    residuals_QR(k) = norm(d);
end


% Obtain the column for which the residual is minimal.

switch SETTINGS.PLOT_GRAPHS
    case 'y'
        [~,opt_col] = min(log10(residuals_QR));
        figure('name','Optimal Column Calculation')
        plot(log10(residuals_QR),'-s');
        hold on
        title('Residuals from removing each column c_{t_{1},t_{2},j} of S_{t_{1},t_{2}}')
        xlabel('k: Index of subresultant colum removednfrom S_{t_{1},t_{2}}')
        ylabel('log_{10} Residual')
        hold off
    case 'n'
    otherwise
        error('err')
end

%% Print the optimal column and minimal residual
% fprintf('\n')
% fprintf('Optimal column for removal is given by %i \n',opt_col)
% fprintf('Minimal Residual %0.5e',min(residuals_QR));
% fprintf('\n')


end