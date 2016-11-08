function [idx_col] = GetOptimalColumn_Total(Sk)
% Get index of the optimal column of the Sylvester Matrix S_{k} to be 
% removed.
%
% % Inputs.
% 
% Sk : Sylvester subresultant matrix S_{k}(f,g)
%
% % Outputs
%
% idx_col : Index of column to be removed


% From the given subresultant find the optimal column for removal.
[~,nCols_T1] = size(Sk);

% Take the QR decomposition of the Subresultant
% [Qk,Rk] = qr(Sk);


vResiduals = zeros(nCols_T1,1);

for i = 1 : 1 : nCols_T1

%     Sk_temp = Sk;
%     % Rem
%     ck = Sk_temp(:,k);
%     [Q,~] = qrdelete(Qk,Rk,k);
%     cd = Q'*ck;
%     d = cd(nCols_T1+1:end,:);
%     residuals_QR(k) = norm(d);
    
    % Get the matrix Ak, which is S_{k} with column c_{k} removed
    Ak = Sk;
    Ak(:,i) = [];
    
    % Get the column c_{k}
    ck = Sk(:,i);
    
    % Solve A*x = b
    x_ls = SolveAx_b(Ak,ck);
    
    % Get the residual from this solution
    vResiduals(i) = norm(ck - (Ak*x_ls));    
    
end

% Obtain the column for which the residual is minimal.
[~,idx_col] = min(log10(vResiduals));

% Obtain the column for which the residual is minimal.
global SETTINGS
switch SETTINGS.PLOT_GRAPHS
    case 'y'
        figure_title = sprintf([mfilename ' : ' 'Optimal Column Calculation' ]);
        [~,idx_col] = min(log10(vResiduals));
        figure('name',figure_title);
        plot(log10(vResiduals),'-s');
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