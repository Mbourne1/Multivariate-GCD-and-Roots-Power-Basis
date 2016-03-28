function C1 = BuildC1(uxy_matrix,t1,t2,m1,m2)
% BuildC1(uxy_matrix,t1,t2,m1,m2)
%
% Build the Matrix C_{1} consisting of coefficients of u(x,y). Used in the
% approximate polynomial factorisation (APF).
% 
% [ C(u) ; C(v)] d = [f;g]
%
%
% %                         Inputs
%
%
% uxy_matrix :  Input quotient polynomial uxy in matrix form. Excluding
%               binomial coefficients. Excluding thetas.
%
% t1 :  Degree of GCD with respect to x
%
% t2 :  Degree of GCD with respect to y
%
% m1 :  Degree of polynomial f with respect to x
%
% m2 :  Degree of polynomial f with respect to y
%


% get size of uxy_matrix
[m1_t1,m2_t2] = GetDegree(uxy_matrix);


% initialise a zero matrix to pad uxy
zero_matrix = zeros(m1+1,m2+1);

nRowsC1 = (m1+1) * (m2+1);
nColsC1 = (t1+1) * (t2+1);

C1 = zeros(nRowsC1,nColsC1);

num_diags = (t1 + 1) + (t2 + 1) + 1;

count = 1;

% for every diagonal of the matrix dxy_mtrx.
for tot = 0:1:num_diags
    
    for i = tot:-1:0
        j = tot-i;

        % if i is within the number of rows of dxy_mtrx
        % if j is within the number of cols of dxy_mtrx
        if (i <= t1)  && (j <= t2)
            
            uxy_matrix_padded = zero_matrix;
            uxy_matrix_padded((i+1):(m1_t1)+(i+1), (j+1):(m2_t2)+(j+1)) = ...
                uxy_matrix;
            
                        
            temp_vec = GetAsVector(uxy_matrix_padded);
            
            C1(:,count) = temp_vec;
            count = count + 1;
        end
    end
    
end


end