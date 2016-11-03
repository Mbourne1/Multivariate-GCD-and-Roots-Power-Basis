function [opt_col] = GetOptimalColumn_Total(Sk)


% From the given subresultant find the optimal column for removal.
[~,nCols_T1] = size(Sk);

% Take the QR decomposition of the Subresultant
[Qk,Rk] = qr(Sk);


residuals_QR = zeros(nCols_T1,1);
for k=1:1:nCols_T1
    Sk_temp = Sk;
    % Rem
    ck = Sk_temp(:,k);
    [Q,~] = qrdelete(Qk,Rk,k);
    cd = Q'*ck;
    d = cd(nCols_T1+1:end,:);
    residuals_QR(k) = norm(d);
end


% Obtain the column for which the residual is minimal.
global SETTINGS
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