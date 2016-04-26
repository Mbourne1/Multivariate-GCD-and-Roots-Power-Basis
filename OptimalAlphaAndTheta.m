function [alpha,theta1,theta2] = OptimalAlphaAndTheta(fxy_matrix, gxy_matrix)
% OptimalAlphaAndTheta(fxy_matrix, gxy_matrix)
%
% Obtain the optimal values of alpha, theta1 and theta2 for the Sylvester 
% matrix of the two polynomials f(x,y) and g(x,y)
% 
% Inputs.
% 
% fxy_matrix : Coefficients of the polynomial f(x,y)
%
% gxy_matrix : Coefficients of the polynomial g(x,y)


% Define vector f
f = [1 -1 0 0 0];

% get the degree of polynomial f and g

[m1,m2] = GetDegree(fxy_matrix);
[n1,n2] = GetDegree(gxy_matrix);

% Assemble the four submatrices of Matrix A
PartOne = zeros((m1+1)*(m2+1),5);
count = 1;
for i1 = 0:1:m1
    for i2 = 0:1:m2
        % Replace all of the 'count-th row'
        PartOne(count,:) = [1 0 -i1 -i2 0];
        count = count + 1 ;
    end
end

PartTwo = zeros((n1+1)*(n2+1),5);
count = 1;
for i1 = 0:1:n1
    for i2 = 0:1:n2
        % Replace all of the 'count-th row'
        PartTwo(count,:) = [1 0 -i1 -i2 -1];
        count = count + 1;
    end
end


PartThree = zeros((m1+1)*(m2+1),5);
count = 1;
for i1 = 0:1:m1
    for i2 = 0:1:m2
        % Replace all of the 'count-th row'
        PartThree(count,:) = [0 -1 i1 i2 0];
        count = count + 1;
    end
end

PartFour = zeros((n1+1)*(n2+1),5);
count = 1;
for i1 = 0:1:n1
    for i2 = 0:1:n2
        % Replace all of the 'count-th row'
        PartFour(count,:) = [0 -1 i1 i2 1];
        count = count + 1;
    end
end


% Now build the vector b = [\lambda ; \mu ;  \rho ;  \tau]
lambda_vec = zeros((m1+1)*(m2+1),1);
count = 1;
for i1 = 0:1:m1
    for i2 = 0:1:m2
        lambda_vec(count) = fxy_matrix(i1+1,i2+1);
        count = count + 1;
    end
end

% lambda_vec2 = GetAsVector(fxy_matrix);


mu_vec = zeros((n1+1)*(n2+1),1);
count = 1;
for i1 = 0:1:n1
    for i2 = 0:1:n2
        mu_vec(count) = gxy_matrix(i1+1,i2+1);
        count = count + 1;
    end
end

% mu_vec = GetAsVector(gxy_matrix);


rho_vec = zeros((m1+1)*(m2+1),1);
count = 1;
for i1 = 0:1:m1
    for i2 = 0:1:m2
        rho_vec(count) = fxy_matrix(i1+1,i2+1);
        count = count + 1;
    end
end

% rho_vec = GetAsVector(fxy_matrix);

tau_vec = zeros((n1+1)*(n2+1),1);
count = 1;
for i1 = 0:1:n1
    for i2 = 0:1:n2
        tau_vec(count) = gxy_matrix(i1+1,i2+1);
        count = count + 1;
    end
end

% tau_vec = GetAsVector(gxy_matrix)



% Get the index of the zero rows in lambda
index_zero_lambda = find(lambda_vec==0);
% Remove the corresponding zeros from lambda
lambda_vec(index_zero_lambda,:) = [];
% Remove the corresponding rows from PartOne Matrix
PartOne(index_zero_lambda,:) = [];


% Get the index of the zero rows in lambda
index_zero_mu = find(mu_vec==0);
% Remove the corresponding zeros from lambda
mu_vec(index_zero_mu,:) = [];
% Remove the corresponding rows from PartOne Matrix
PartTwo(index_zero_mu,:) = [];


% Get the index of the zero rows in rho
index_zero_rho = find(rho_vec==0);
% Remove the corresponding zeros from rho
rho_vec(index_zero_rho,:) = [];
% Remove the corresponding rows from PartOne Matrix
PartThree(index_zero_rho,:) = [];


% Get the index of the zero rows in tau
index_zero_tau = find(tau_vec==0);
% Remove the corresponding zeros from tau
tau_vec(index_zero_tau,:) = [];
% Remove the corresponding rows from PartOne Matrix
PartFour(index_zero_tau,:) = [];

A =-[PartOne; PartTwo; PartThree; PartFour];


b = -[log10(abs(lambda_vec)); log10(abs(mu_vec)); -log10(abs(rho_vec));-log10(abs(tau_vec))];


x = linprog(f,A,b);

try
    theta1 = 10^x(3);
    theta2 = 10^x(4);
    alpha = 10^x(5);
    %fprintf('Optimal theta 1 and theta 2 given by: \n  theta_{1}: %0.5e \n  theta_{2}: %0.5e',theta1,theta2)
catch
    fprintf('Failed to optimize\n')
    theta1 = 1;
    theta2 = 1;
    alpha = 1; 

end


end