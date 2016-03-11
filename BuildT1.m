function T1 = BuildT1(fxy_matrix,n1_k1,n2_k2)
% Given the polynomial f(x,y) build the partition T1 of the Sylvester
% matrix S(f,g) = [T(f) T(g)]. 
%
%           Inputs.
%
% fxy_matrix :
%
% n1 :
%
% n2 :
%
% k1 :
%
% k2 :

%% Get structure of f(x,y)
% Get size of fxy
[rows,cols] = size(fxy_matrix);

% Get degree of f with respect to x
m1 = rows - 1;
% Get degree of f with respect to y
m2 = cols - 1;

%%
% Get number of multiplications with respect to x
mult_wrt_x = n1_k1 + 1;

% Get number of multiplications with respect to y
mult_wrt_y = n2_k2 + 1;

% Initalise a zero matrix
zero_matrix = zeros(m1 + n1_k1 + 1, m2 + n2_k2 + 1);

T1 = [];

% Get number of diagonals in the matrix v(x,y)
num_diags = (n1_k1+1) + (n2_k2+1) - 1;

% for each diagonal 
for tot = 0 : 1 : num_diags
    for i = tot:-1:0
        
        % Initialise a temporary zero matrix
        temp_mat = zero_matrix;
        
        j = tot - i;
        
        if i <= (n1_k1) && j <= (n2_k2)
            % Wrap the entries of fxy_matrix_padded 1 downward
            % wrap the entires of fxy_matrix 1 place to the right
                       
            temp_mat((i+1):(rows+i),(j+1):(cols+j)) = fxy_matrix;
            
            
            % Produce temporary vector from the coefficients       
            temp_vec = GetAsVector(temp_mat);
            
            % Append the temporary vector to the matrix T1
            T1 = [T1 temp_vec];
        end
    end
    
end
end