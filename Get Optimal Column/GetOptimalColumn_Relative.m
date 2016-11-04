function [idx_col] = GetOptimalColumn_Relative(Sk)
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
[~,nCols_T1] = size(Sk);

% Take the QR decomposition of the Subresultant
[Qk,Rk] = qr(Sk);

n = nCols_T1;


residuals_QR = zeros(nCols_T1,1);
for k=1:1:nCols_T1
    Sk_temp = Sk;
    % Rem
    ck = Sk_temp(:,k);
    [Q,~] = qrdelete(Qk,Rk,k);
    cd = Q'*ck;
    d = cd(n+1:end,:);
    residuals_QR(k) = norm(d);
end


% Obtain the column for which the residual is minimal.
[~,idx_col] = min(log10(residuals_QR));

try
switch SETTINGS.PLOT_GRAPHS
    case 'y'       
        figure('name','Optimal Column Calculation')
        plot(log10(residuals_QR),'-s');
        hold on
        title('Residuals from removing each column c_{t_{1},t_{2},j} of S_{t_{1},t_{2}}')
        xlabel('k: Index of subresultant colum removednfrom S_{t_{1},t_{2}}')
        ylabel('log_{10} Residual')
        hold off
    case 'n'
    otherwise
        error('Error: plot_graphs is either y or n')
end
catch err
    fprintf(err.message)
end
%% Print the optimal column and minimal residual
% fprintf('\n')
% fprintf('Optimal column for removal is given by %i \n',opt_col)
% fprintf('Minimal Residual %0.5e',min(residuals_QR));
% fprintf('\n')


end