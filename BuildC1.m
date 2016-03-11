function C1 = BuildC1(uxy_matrix,t1,t2,m1,m2)
% Build the Matrix C_{1} consisting of coefficients of u(x,y). Used in the
% approximate polynomial factorisation (APF).
% 
% [ C(u) ; C(v)] d = [f;g]
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% get size of uxy_matrix
[rows,cols] = size(uxy_matrix);

% get m1-t1, degree of uxy with respect to x
m1_t1 = rows - 1;

% Get m2-t2, degree of uxy with respect to y
m2_t2 = cols - 1;

% initialise a zero matrix to pad uxy
zero_matrix = zeros(m1+1,m2+1);

%
C1 = [];

num_diags = (t1 + 1) + (t2 + 1) + 1;

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
            
            C1 = [C1 temp_vec];
        end
    end
    
end