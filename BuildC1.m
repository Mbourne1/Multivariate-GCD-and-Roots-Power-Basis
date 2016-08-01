function C1 = BuildC1(uxy_matrix,t1,t2)
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



% Get the degree of u(x,y)
[m1_t1,m2_t2] = GetDegree(uxy_matrix);

% Get m1 and m2
m1 = m1_t1 + t1;
m2 = m2_t2 + t2;

% initialise a zero matrix to pad uxy
zero_matrix = zeros(m1+1,m2+1);

nCoefficients_fxy = (m1+1) * (m2+1);
nCoefficients_dxy = (t1+1) * (t2+1);


nRowsC1 = nCoefficients_fxy;
nColsC1 = nCoefficients_dxy;

C1 = zeros(nRowsC1,nColsC1);

% Get number of diagonals in the matrix of coefficients of d(x,y)
nDiags_dxy = (t1 + 1) + (t2 + 1) + 1;

count = 1;

% For every diagonal of the matrix of coefficients of d(x,y).
for tot = 0:1:nDiags_dxy
    
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