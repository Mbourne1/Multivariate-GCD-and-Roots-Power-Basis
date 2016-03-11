
function T1 = BuildT1_totaldegree(fxy_matrix,m,n,k)

% get size of f(x,y)
[rows,cols] = size(fxy_matrix);

% Initalise a zero matrix
zero_matrix = zeros(m+n-k+1, m+n-k+1);

T1 = [];

% Get number of diagonals in the matrix v(x,y)
num_diags = (n-k+1);

% for each diagonal in M(v)
for tot = 0 : 1 : num_diags -1
    for i = tot:-1:0
        
        % Initialise a temporary zero matrix
        temp_mat = zero_matrix;
        
        j = tot - i;
        
        if i <= (n-k) && j <= (n-k)
            % Wrap the entries of fxy_matrix_padded 1 downward
            % wrap the entires of fxy_matrix 1 place to the right
                       
            temp_mat((i+1):(rows+i),(j+1):(cols+j)) = fxy_matrix;
            
            % Produce temporary vector from the coefficients       
            temp_vec = GetAsVector(temp_mat);
            
            % remove the last nchoosek(n-k+1,2) entries from temp vec
            % only keep the nchoosek(n-k+2,2) non-zero values
            val = nchoosek(m+n-k+2,2);
            temp_vec2 = temp_vec(1:nchoosek(m+n-k+2,2));
            
            T1 = [T1 temp_vec2];
        end
    end
    
end
end