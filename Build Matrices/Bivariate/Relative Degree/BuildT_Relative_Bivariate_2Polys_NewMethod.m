function T = BuildT_Relative_Bivariate_2Polys_NewMethod(fxy, gxy, k1, k2)
% Build the Sylvester matrix 
%
% % Inputs
%
% fxy : Coefficients of polynomial f(x,y)
% 
% gxy : Coefficients of polynomial g(x,y)
%
% k1
%
% k2

% Get degree of f(x,y)
[m1,m2] = GetDegree_Bivariate(fxy);

% Get degree of g(x,y)
[n1,n2] = GetDegree_Bivariate(gxy);

T1 = BuildT1_Relative_Bivariate_2Polys_NewMethod(fxy, n1-k1, n2-k2);
T2 = BuildT1_Relative_Bivariate_2Polys_NewMethod(gxy, m1-k1, m2-k2);

T = [T1 T2];


    
end

function T1 = BuildT1_Relative_Bivariate_2Polys_NewMethod(fxy,n1_k1,n2_k2)

[m1, m2] = GetDegree_Bivariate(fxy);

% Get array of polynomials f_{i}(x)
arr_fx = cell(m2+1,1);
for i = 0:1:m2
   arr_fx{i+1} = fxy(:,i+1); 
end


% Get array of partitions of the Sylvester matrix S_{k1,k2}(f(x,y),g(x,y))
arr_T1 = cell(m2 + n2_k2 + 1, n2_k2 + 1);
[arr_T1{:}] = deal(zeros(m1+n1_k1 + 1, n1_k1+1));

for j = 0:1:n2_k2
    for i = j : 1 : j + m2
        
        arr_T1{i+1,j+1} = BuildT1_Univariate(arr_fx{i-j+1}, n1_k1);
        
    end
end

T1 = cell2mat(arr_T1);

end

