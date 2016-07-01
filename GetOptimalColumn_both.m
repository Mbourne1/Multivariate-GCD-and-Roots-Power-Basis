function [opt_col] = GetOptimalColumn_both(fxy_matrix,gxy_matrix,m,n,t,t1,t2)
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

% Get degree structure of polynomial f(x,y)
[m1,m2] = GetDegree(fxy_matrix);

% Get degree structure of polynomial g(x,y)
[n1,n2] = GetDegree(gxy_matrix);

% % Build the partitions of the Sylvester matrix S_{t}

% Build the first partition containing coefficients of fxy
T1 = BuildT1(fxy_matrix,n1-t1,n2-t2);

% Build the second partition containing coefficients of gxy
T2 = BuildT1(gxy_matrix,m1-t1,m2-t2);

% Remove the columns of T1 corresponding to the zeros in v(x,y)
nNonZeros_vxy = GetNumNonZeros(n1-t1,n2-t2,n-t);
T1 = T1(:,1:nNonZeros_vxy);

% Remove the columns of T2 corresponding to the zeros in u(x,y)
nNonZeros_uxy = GetNumNonZeros(m1-t1,m2-t2,m-t);
T2 = T2(:,1:nNonZeros_uxy);

% %
% %
% Remove the rows of T1
nNonZeros_fv = GetNumNonZeros(m1+n1-t1,m2+n2-t2,m+n-t);
nNonZeros_gu = GetNumNonZeros(n1+m1-t1,n2+m2-t2,n+m-t);

T1 = T1(1:nNonZeros_fv,:);
T2 = T2(1:nNonZeros_gu,:);

% Concatenate the two partitions
S_t1t2 = [T1 T2];

% From the given subresultant find the optimal column for removal.
[~,nCols_S] = size(S_t1t2);



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
    Ak = S_t1t2;
    Ak(:,i) = [];
    
    % Get the column c_{k}
    ck = S_t1t2(:,i);
    
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