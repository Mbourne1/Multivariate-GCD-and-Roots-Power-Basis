
function T1 = BuildT1_Total_Bivariate_NewMethod(fxy, m, n_k)


% Get array of polynomails f_{i}(x,y) where each f_{i}(x,y) is of degree i.
for i = 0:1:m
   
   % Initialise polynomial f_{i}(x,y)
   temp_poly = zeros(i+1,1); 
    
   % Get coefficients of polynomial f_{i}(x,y) of degree i
   for j = 0:1:i
       temp_poly(j+1) = fxy(i-j+1,j+1);
   end
   
   arr_fxy{i+1} = temp_poly;
   
end


% Build the convolution matrix of size n-k
arr = cell(m+n_k + 1, n_k+1);

for i = 0:1:m + n_k
    for j = 0:1: n_k
        
        if i-j>= 0 && i-j < m
            arr{i+1,j+1} = BuildT1_Univariate(arr_fxy{i-j+1},j);
        else
            arr{i+1,j+1} = zeros(i+1, j+1);
        end
    end
    
end

T1 = cell2mat(arr);

end