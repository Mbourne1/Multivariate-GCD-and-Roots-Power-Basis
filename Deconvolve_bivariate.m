function [hxy_matrix] = Deconvolve_bivariate(fxy_matrix,gxy_matrix)
% return the matrix of coefficients of the polynomial h, where h = f/g

%%

% Get degrees of polynomial f
[r,c] = size(fxy_matrix);
m1 = r - 1;
m2 = c - 1;

% Get the degrees of polynomial g
[r,c] = size(gxy_matrix);
n1 = r - 1;
n2 = c - 1;

% multiply gxy by the (m1-n1+1)*(m2-n2+1)

% build the matrix which will contain the multiplied entries
zero_matrix = zeros(m1+1,m2+1);

% the number of diagonals of the multiplication matrix h
num_diags = (m1 - n1 + 1) + (m2 - n2 + 1) + 1;

C1 = [];
% for every diagonal of the matrix dxy_mtrx.
for tot = 0:1:num_diags
    
    for i = tot:-1:0
        j = tot-i;

        % if i is within the number of rows of hxy_mtrx
        % if j is within the number of cols of hxy_mtrx
        if (i <= m1 - n1)  && (j <= m2 - n2)
            
            gxy_matrix_padded = zero_matrix;
            gxy_matrix_padded((i+1):(n1)+(i+1), (j+1):(n2)+(j+1)) = ...
                gxy_matrix;
            
                        
            temp_vec = GetAsVector(gxy_matrix_padded);
            
            C1 = [C1 temp_vec];
        end
    end
    
end

%% Get the polynomial f in vector form
f = GetAsVector(fxy_matrix);

h = pinv(C1) * f;

hxy_matrix = GetAsMatrix(h,m1-n1,m2-n2);

hxy_matrix;

