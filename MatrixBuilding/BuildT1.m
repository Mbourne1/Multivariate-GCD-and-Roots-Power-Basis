function T1 = BuildT1(fxy_matrix,n1_k1,n2_k2)
% Given the polynomial f(x,y) build the partition T(f) of the Sylvester
% matrix S(f,g) = [T(f) T(g)]. 
% T() * v = T(g) * u, so [T(f) T(g)] * [v -u] = 0
%
% BuildT1(fxy_matrix,n1_k1,n2_k2)
%
% Inputs.
%
% fxy_matrix : Coefficient matrix of f(x,y)
%
% n1_k1 : Degree of v(x,y) with respect to x
%
% n2_k2 : Degree of v(x,y) with respect to y
% 
% Outputs.
%
% T1 : Partiton T_{n1-k1,n2-k2}(f) of the Sylvester matrix

% Get degree of polynomial f(x,y) with respect to x and y
[m1,m2] = GetDegree(fxy_matrix);

[nRows_f,nCols_f] = size(fxy_matrix);

% The coefficients of f(x,y) must be multiplied by the basis elements of
% v(x,y) to form the columns of S_{k_{1},k_{2}}(f,g). 

% Initalise a zero matrix
zero_matrix = zeros(m1 + n1_k1 + 1, m2 + n2_k2 + 1);

nRows_T1 = (m1 + n1_k1 + 1) * (m2 + n2_k2 + 1);
nCols_T1 = (n1_k1 + 1) * (n2_k2 + 1);

T1 = zeros(nRows_T1,nCols_T1);

% Get number of diagonals in the matrix v(x,y)
nDiags_v = (n1_k1+1) + (n2_k2+1) - 1;

count = 1;
% for each diagonal in v(x,y)
for tot = 0 : 1 : nDiags_v
    
    for i = tot:-1:0
        
        % Initialise a temporary zero matrix
        temp_mat = zero_matrix;
        
        j = tot - i;
        
        if i <= (n1_k1) && j <= (n2_k2)
            % Wrap the entries of fxy_matrix_padded 1 downward
            % wrap the entires of fxy_matrix 1 place to the right
                       
            temp_mat((i+1):(nRows_f+i),(j+1):(nCols_f+j)) = fxy_matrix;
            
            
            % Produce temporary vector from the coefficients       
            temp_vec = GetAsVector(temp_mat);
            
            % Set the kth column 
            
            T1(:,count) = temp_vec;
            count = count + 1;
        end
    end
    
end


end


