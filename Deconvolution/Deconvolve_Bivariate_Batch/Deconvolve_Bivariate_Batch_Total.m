function arr_hxy = Deconvolve_Bivariate_Batch_Total(arr_fxy,vDegt_arr_fxy)
% Perform the batch deconvolution C(f1,f2,..,fn-1) h = [f0;f1;,...,;fn]
%
% Inputs.
%
% arr_fxy : Array of polynomials f(x,y)
%
% vDeg_fxy : Vector containing the total degree of the polynomials f{i}

% %
% %
% Form the left hand side matrix

% Get number of polynomials in the array arr_fxy
nPolys_arr_fxy = size(arr_fxy,1);

vDegt_arr_hxy = abs(diff(vDegt_arr_fxy));


% For each of the polynomials excluding the first.
for i = 2:1:nPolys_arr_fxy 
    
     
    % Build the matrix T(f(x,y))
    T1 = BuildT1_Total(arr_fxy{i,1},vDegt_arr_fxy(i,1),vDegt_arr_hxy(i-1,1));
    
    % % Strip the columns corresponding to zeros in f{i-1}
    arr_T1{i-1} = T1;
end

C = blkdiag(arr_T1{:});

% %
% %
% Form the Right hand side vector

for i = 1:1:nPolys_arr_fxy - 1
    
    % Temporarily label the ith entry of the array as fxy
    v_fxy = GetAsVector(arr_fxy{i});
    
    % Remove the corresponding zeros
    
    % Get number of zeros
    m = vDegt_arr_fxy(i);
    nNonZeros_fxy = nchoosek(m+2,2);
    
    % Remove zeros from vector
    v_fxy = v_fxy(1:nNonZeros_fxy);
    
    arr_rhs{i} = v_fxy;
    
end

vRHS = cell2mat(arr_rhs');

% %
% %
% %
x_ls = SolveAx_b(C,vRHS);

% Split solution vector into polynomials h{i}.
nPolys_hxy = nPolys_arr_fxy -1 ;

arr_hxy = cell(nPolys_hxy,1);

for i = 1:nPolys_hxy
    
    
    n = vDegt_arr_hxy(i);
    nNonZeros = nchoosek(n+2,2);
    
    % Get vector of zeros 
    temp_vec = zeros((n+1)*(n+1),1);
    temp_vec(1:nNonZeros) = x_ls(1:nNonZeros);
    
    mat = GetAsMatrix(temp_vec,n,n);
    
    arr_hxy{i} = mat;
    x_ls(1:nNonZeros) = [];
    
end
end