function [opt_col] = GetOptimalColumn_Respective(fxy_matrix,gxy_matrix,t1,t2)
% Given two polynomials f(x,y) and g(x,y), and the degree structure of
% their GCD, pick the optimal column for removal from the Sylvester matrix
% S_{t_{1},t_{2}}(f,g).

global SETTINGS

% Get degree structure of polynomial f(x,y)
[m1,m2] = GetDegree(fxy_matrix);

% Get degree structure of polynomial g(x,y)
[n1,n2] = GetDegree(gxy_matrix);

% % Build the partitions of the Sylvester matrix S_{t}

% Build the first partition containing coefficients of fxy
T1 = BuildT1(fxy_matrix,n1-t1,n2-t2);

% Build the second partition containing coefficients of gxy
T2 = BuildT1(gxy_matrix,m1-t1,m2-t2);

% Concatenate the two partitions
St = [T1 T2];

% From the given subresultant find the optimal column for removal.
[~,cols_T1] = size(St);

% Take the QR decomposition of the Subresultant
[Qk,Rk] = qr(St);


n = n1+n2;

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
[~,opt_col] = min(log10(residuals_QR));

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