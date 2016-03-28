function C1 = BuildC1_total(uxy_matrix,m,t)
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
% t1 : Total Degree of GCD 
%
% m :  Total Degree of polynomial f
%
% n :  Total Degree of polynomial f.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% initialise a zero matrix to pad uxy
zero_matrix = zeros(m+1,m+1);

%
C1 = [];

% Get number of diagonals in the matrix M(d(x,y))
num_diags_d = (t + 1); 

% for every diagonal of the matrix dxy_mtrx.
for tot = 0:1:num_diags_d
    
    for i = tot:-1:0
        j = tot-i;

        % if i is within the number of rows of dxy_mtrx
        % if j is within the number of cols of dxy_mtrx
        if (i <= t)  && (j <= t)
            
            uxy_matrix_padded = zero_matrix;
            uxy_matrix_padded((i+1):(m-t)+(i+1), (j+1):(m-t)+(j+1)) = ...
                uxy_matrix;
            
                        
            temp_vec = getAsVector(uxy_matrix_padded);
            % Remove the final nchoosek(m+2-1,2)
            temp_vec = temp_vec(1:nchoosek(m+2,2));
           
            C1 = [C1 temp_vec];
        end
    end
    
end

end