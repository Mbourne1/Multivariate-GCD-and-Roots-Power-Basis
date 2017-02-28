function [idx_col] = GetOptimalColumn_Relative(Sk1k2)
% Given two polynomials f(x,y) and g(x,y), and the degree structure of
% their GCD, pick the optimal column for removal from the Sylvester matrix
% S_{t_{1},t_{2}}(f,g).
%
% % Inputs.
%
% Sk : Sylvester Subresultant matrix S_{k_{1},k_{2}}
%
% % Outputs.
%
% idx_col : Index of optimal column for removal

global SETTINGS


% From the given subresultant find the optimal column for removal.
[~,nCols_Sk] = size(Sk1k2);

% Take the QR decomposition of the Sylvester subresultant
% % [Qk,Rk] = qr(Sk);



vResiduals= zeros(nCols_Sk,1);

for i = 1 : 1 : nCols_Sk
    
    %     Sk_temp = Sk;
    %     ck = Sk_temp(:,k);
    %     [Q,~] = qrdelete(Qk,Rk,k);
    %     cd = Q'*ck;
    %     d = cd(n+1:end,:);
    %     vResiduals(k) = norm(d);
    % Get the matrix Ak, which is S_{k} with column c_{k} removed
    
    Ak = Sk1k2;
    Ak(:,i) = [];
    
    % Get the column c_{k}
    ck = Sk1k2(:,i);
    
    % Solve A*x = b
    x_ls = SolveAx_b(Ak,ck);
    
    % Get the residual from this solution
    vResiduals(i) = norm(ck - (Ak*x_ls));
end


% Obtain the column for which the residual is minimal.
[~,idx_col] = min(log10(vResiduals));

% Obtain the column for which the residual is minimal.
if(SETTINGS.PLOT_GRAPHS)
    figure_title = sprintf([mfilename ' : ' 'Optimal Column Calculation' ]);
    [~,idx_col] = min(log10(vResiduals));
    figure('name',figure_title);
    plot(log10(vResiduals),'-s');
    hold on
    title('Residuals from removing each column c_{t_{1},t_{2},j} of S_{t_{1},t_{2}}')
    xlabel('k: Index of subresultant colum removednfrom S_{t_{1},t_{2}}')
    ylabel('log_{10} Residual')
    hold off
end

end