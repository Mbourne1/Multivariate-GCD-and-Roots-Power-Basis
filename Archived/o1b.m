function [uxy_mtrx_clc, vxy_mtrx_clc, dxy_mtrx_clc] = o1(fxy_matrix_working,gxy_matrix_working,...
    m,n)
global bool_preproc

% fxy_matrix_before_proc = fxy_matrix_exact;
% gxy_matrix_before_proc = gxy_matrix_exact;

%% Preprocessing
switch bool_preproc
    case 'y'
        % Try to calculate geometric mean of entries
        try
            GM_f = geomean(abs(nonzeros(fxy_matrix_working)));
            GM_g = geomean(abs(nonzeros(gxy_matrix_working)));
        catch
            GM_f = 1;
            GM_g = 1;
        end
        
        % Normalise coefficients of f and g by dividing by geometric
        % mean
        fxy_matrix_working = fxy_matrix_working./ GM_f;
        gxy_matrix_working = gxy_matrix_working./ GM_g;
        
        
        [opt_alpha, opt_theta_1,opt_theta_2] = OptimalAlphaAndTheta(fxy_matrix_working,gxy_matrix_working);
                
        
        %[opt_theta_1,opt_theta_2] = OptimalTheta(fxy_matrix_working,gxy_matrix_working);
        
        
        fprintf('Optimal theta_{1}  :  %0.5e \n',opt_theta_1)
        fprintf('Optimal theta_{2}  :  %0.5e \n',opt_theta_2)
        fprintf('Optimal alpha :   %0.5e \n', opt_alpha)
        
        % multiply each of the coefficients by the optimal theta
        
        % multiply fxy by optimal theta 1
        % each row of the matrix is multiplied by theta_{1}^{i} where i is
        [r,c] = size(fxy_matrix_working);
        
        for i = 1:1:r
            fxy_matrix_working(i,:) = fxy_matrix_working(i,:) .* (opt_theta_1.^(i-1)) ;
        end
        for j = 1:1:c
            fxy_matrix_working(:,j) = fxy_matrix_working(:,j) .* (opt_theta_2.^(j-1)) ;
        end
        
        
        % each row of the matrix is multiplied by theta_{1}^{i} where i is
        [r,c] = size(gxy_matrix_working);
        
        for i = 1:1:r
            gxy_matrix_working(i,:) = gxy_matrix_working(i,:) .* (opt_theta_1.^(i-1)) ;
        end
        for j = 1:1:c
            gxy_matrix_working(:,j) = gxy_matrix_working(:,j) .* (opt_theta_2.^(j-1)) ;
        end
        
        % multiply gxy by alpha
        gxy_matrix_working = gxy_matrix_working .*opt_alpha;
               
    case 'n'
        opt_theta_1 = 1;
        opt_theta_2 = 1;
        opt_alpha = 1;
end


% Having obtained preprocessed versions, compare the max:min of unprocessed
% polynomials, to max:min preprocessed.



%% Get the Degree of the GCD
[t1,t2] = getGCDDegree(fxy_matrix_working,gxy_matrix_working,m,n,opt_theta_1, opt_theta_2);
t = t1 + t2;


fprintf('----------------------------------------------------------------\n')
fprintf('\n')
fprintf('The Calculated Degree of the GCD is given by \n')
fprintf('Degree of GCD wrt x : t1 = %i\n',t1)
fprintf('Degree of GCD wrt y : t2 = %i\n',t2)
fprintf('\n')
fprintf('----------------------------------------------------------------\n')


% get the degrees of the input polynomials
[r,c] = size(fxy_matrix_working);
m1 = r-1;
m2 = c-1;

%
[r,c] = size(gxy_matrix_working);
n1 = r-1;
n2 = c-1;

%% Get the quotient polynomials uxy and vxy

[uxy_mtrx_clc,vxy_mtrx_clc] = getQuotients(fxy_matrix_working,gxy_matrix_working,m,n,t,t1,t2,opt_theta_1,opt_theta_2);

% strip alpha from the coefficients of uxy
uxy_mtrx_clc = uxy_mtrx_clc ./ opt_alpha;

% Normlise by dividing by the first non-zero coefficient which is given by
% the the coefficient u_{i1,i2}x^{t1}y^{t2}
% where:
%   i_{1} = m_{1}-t_{1}+1;
%   i_{2} = m_{2}-t_{2}+1;

% Normalise


%% Get the GCD dxy

dxy_mtrx_clc = getGCDCoefficeints(fxy_matrix_working,gxy_matrix_working,uxy_mtrx_clc,vxy_mtrx_clc, opt_theta_1, opt_theta_2);

% Strip thetas from uxy_matrix vxy_matrix and dxy_matrix
switch bool_preproc
    case 'y'
        [r,c] = size(uxy_mtrx_clc);
        % divide the first row by theta_{1}^{0}
        for i=1:1:r
            uxy_mtrx_clc(i,:) = uxy_mtrx_clc(i,:) ./ (opt_theta_1^(i-1));
        end
        for i=1:1:c
            uxy_mtrx_clc(:,i) = uxy_mtrx_clc(:,i) ./ (opt_theta_2^(i-1)); 
        end
        
        
        [r,c] = size(vxy_mtrx_clc);
        % divide the first row by theta_{1}^{0}
        for i=1:1:r
            vxy_mtrx_clc(i,:) = vxy_mtrx_clc(i,:) ./ (opt_theta_1^(i-1));
        end
        for i=1:1:c
            vxy_mtrx_clc(:,i) = vxy_mtrx_clc(:,i) ./ (opt_theta_2^(i-1)); 
        end
         
        
        [r,c] = size(dxy_mtrx_clc);
        % divide the first row by theta_{1}^{0}
        for i=1:1:r
            dxy_mtrx_clc(i,:) = dxy_mtrx_clc(i,:) ./ (opt_theta_1^(i-1));
        end
        for i=1:1:c
            dxy_mtrx_clc(:,i) = dxy_mtrx_clc(:,i) ./ (opt_theta_2^(i-1)); 
        end
        
end


fprintf('----------------------------------------------------------------\n')
fprintf('\n')
fprintf('Calculated values:\n')
fprintf('uxy_matrix\n')
uxy_mtrx_clc ./ uxy_mtrx_clc(1,1)


fprintf('\n')
fprintf('----------------------------------------------------------------\n')
fprintf('\n')
fprintf('vxy_matrix\n')
vxy_mtrx_clc ./ vxy_mtrx_clc(1,1)

fprintf('\n')
fprintf('----------------------------------------------------------------\n')

fprintf('dxy_matrix\n')
dxy_mtrx_clc./dxy_mtrx_clc(1,1)


%% Testing the result
%multiply u by d to see if we obtain coefficients of f

% Build the cauchy matrix C(u)
C_u = BuildC1(uxy_mtrx_clc,t1,t2,m1,m2,1,1);

% Build the cauchy matrix C(v)
C_v = BuildC1(vxy_mtrx_clc,t1,t2,n1,n2,1,1);

% Build the vector d from the matrix d
count = 1;
for tot = 0:1:(t1 + t2+1)
    for ihat = 0:1:tot
        jhat = tot - ihat;
        if ihat <= t1 && jhat <= t2
            dxy_vec_calc(count,1) = dxy_mtrx_clc(ihat+1,jhat+1);
            count = count + 1;
        end
    end
end
