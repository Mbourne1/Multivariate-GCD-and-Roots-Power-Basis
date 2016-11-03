function [opt_col] = GetOptimalColumn_Both(S_kk1k2)
% Given two polynomials f(x,y) and g(x,y), and the degree structure of
% their GCD t (t1,t2), pick the optimal column for removal from the 
% Sylvester matrix S_{t_{1},t_{2}}(f,g), where some columns are removed
% from each partition of S_{t1,t2}(f,g) corresponding to the zeros in
% u(x,y) and v(x,y).
%
% Inputs.
%
% fxy_maatrix : Coefficients of f(x,y)
%
% gxy_matrix : Coefficients of g(x,y)
%
% t : total degree of d(x,y)
%
% t1 : degree of d(x,y) with respect to x
%
% t2 : degree of d(x,y) with respect to y
%
% Outputs.
%
% opt_col : Index of optimal column for removal 

global SETTINGS


% From the given subresultant find the optimal column for removal.
[~,nCols_S] = size(S_kk1k2);



% Initialise a vector to store residuals by QR decomposition
residuals_QR = zeros(nCols_S,1);
vResiduals = zeros(nCols_S,1);


% % Take the QR decomposition of the Subresultant
% [Qk,Rk] = qr(S_t1t2);
% 
% 
% n = n1+n2;
% for k=1:1:nCols_S
%     
%     % Initialise the temporary Sylvester matrix
%     S_temp = S_t1t2;
%     
%     % Remove column c_{k} from S_{t1t2}
%     ck = S_temp(:,k);
%     
%     % Update QR decomposition of S_{t1t2} with column removed
%     [Q,~] = qrdelete(Qk,Rk,k);
%     
%     cd = Q'*ck;
%     d = cd(n+1:end,:);
%     residuals_QR(k) = norm(d);
% end


for i = 1:1:nCols_S
    
    % Get the matrix Ak, which is S_{k} with column c_{k} removed
    Ak = S_kk1k2;
    Ak(:,i) = [];
    
    % Get the column c_{k}
    ck = S_kk1k2(:,i);
    
    % Solve A*x = b
    x_ls = SolveAx_b(Ak,ck);
    
    % Get the residual from this solution
    vResiduals(i) = norm(ck - (Ak*x_ls));
    
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