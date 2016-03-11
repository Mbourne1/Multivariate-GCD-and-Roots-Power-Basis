
function T1 = BuildT1new(fxy_matrix,n1,n2,k1,k2)

% get size of fxy
[rows,cols] = size(fxy_matrix);

% Get degree of f with respect to x
m1 = rows - 1;

% Get degree of f with respect to y
m2 = cols - 1;

% Get number of multiplications with respect to x
mult_wrt_x = n1 - k1 + 1;

% Get number of multiplications with respect to y
mult_wrt_y = n2 - k2 + 1;

% Initalise a zero matrix
zero_matrix = zeros(m1 + n1 - k1 + 1, m2 +n2 - k2 + 1);

T1 = [];
% for each diagonal
for tot = 0:1:(n2-k2+1)+(n1-k1+1) -1;
%fprintf('The diagonal, index : %i \n' , tot)    
    for j = 0:1:tot
        
        % Initialise a temporary zero matrix
        temp_mat = zero_matrix;
        
        i = tot - j;
        
        if i <= n1-k1 && j <= n2-k2
            % Wrap the entries of fxy_matrix_padded 1 downward
            % wrap the entires of fxy_matrix 1 place to the right
                       
            temp_mat((i+1):(rows+i),(j+1):(cols+j)) = fxy_matrix;
            
            % Produce temporary vector from the coefficients
            
            [r,c] = size(temp_mat);
            temp_vec = [];
            for tot2 = 0 : 1: r+c-2
                for j2 = 0:1:tot2
                    
                    i2 = tot2 - j2;
                    
                    if (i2<r) && (j2 < c)
                        temp_vec = [temp_vec ; temp_mat(i2+1,j2+1)];
                    end
                end
            end
            
            T1 = [T1 temp_vec];
        end
    end
    
end
end