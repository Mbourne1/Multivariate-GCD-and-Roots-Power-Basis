
function T1 = BuildT1(fxy_matrix,n1,n2,k1,k2,th1,th2)

% get size of fxy
[rows,cols] = size(fxy_matrix);

% Get degree of f with respect to x
m1 = rows - 1;

% Get degree of f with respect to y
m2 = cols - 1;

% Build the matrix which pre multiplies f_xy
Pre_thetas = zeros(m1+1,1);
for i = 0:1:m1
    Pre_thetas(i+1) = th1^i;
end
Pre_thetas_mtrx = diag(Pre_thetas);

% Build the matrix which post multiplies f_xy
Post_thetas = zeros(m2+1,1);
for i = 0:1:m2
   Post_thetas(i+1) = th2^i; 
end
Post_thetas_mtrx = diag(Post_thetas);

% Create matrix fxy_bi which includes the binomial coefficients
fxy_matrix_new =  fxy_matrix;

% Include thetas and binomial coefficients
fw1w2_matrix_new = Pre_thetas_mtrx * fxy_matrix_new * Post_thetas_mtrx;

% Get number of multiplications with respect to x
mult_wrt_x = n1 - k1 + 1;

% Get number of multiplications with respect to y
mult_wrt_y = n2 - k2 + 1;

% Initalise a zero matrix 
zero_matrix = zeros(m1 + n1 - k1 + 1, m2 +n2 - k2 + 1);
fxy_matrix_padded = zero_matrix;

%input fxy into the first segment of the padded zero matrix 
fxy_matrix_padded(1:m1+1,1:m2+1) = fw1w2_matrix_new;

T1 = [];

% for every diagonal of the matrix vxy. 
for tot0 = 0:1:(n2-k2)+(n1-k1)
    
    for j1 = tot0:-1:0
        j2 = tot0-j2;
        % if j1 is within the number of columns of vxy_mtrx.
        % if j2 is within the number of rows of vxy_mtrx.
        if j1 < n2-k2 + 1 && j2 < n1-k1 +1
            
            %Initialise a zero matrix
            fxy_matrix_bi_padded_new = zero_matrix;
            
            % place the coefficients of fxy_matrix into the zero matrix
            fxy_matrix_bi_padded_new((j2+1):(m1)+(j2+1),(j1+1):(m2)+(j1+1)) = fw1w2_matrix_new;
            
            % convert to a vector
            temp_vec = zeros(m1+n1-k1+1 * m2+n2-k2+1,1);
            count = 1;
            for tot = 0:1:(m1+m2+n1+n2-k1-k2)
                for i2 = tot:-1:0
                    i1 = tot - i2;
                    if i2 <= m2+n2-k2 && i1 <= m1+n1-k1
                        temp_vec(count) = fxy_matrix_bi_padded_new(i1+1,i2+1) ;
                        count = count+1;
                    end
                end
            end
            
            T1 = [T1 temp_vec];
        end
        %
    end
end

T1 = T1;

end