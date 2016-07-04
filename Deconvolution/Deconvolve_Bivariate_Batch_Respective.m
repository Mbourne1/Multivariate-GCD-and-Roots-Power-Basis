function arr_hxy = Deconvolve_Bivariate_Batch_Respective(arr_fxy,vDeg_fxy)
% Perform the batch deconvolution
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
nPolys_fxy = size(arr_fxy,2);

% For each of the polynomials excluding the first.
for i = 2:1:nPolys_fxy
    
    % Get degree of f{i-1} (The corresponding vector in the RHS)
    [n1,n2] = GetDegree(arr_fxy{i-1});
    
    % Temporarily call the ith entry f(x,y)
    fxy = arr_fxy{i};
    
    % Get degree of f(x,y) with respect to x and y
    [m1,m2] = GetDegree(fxy);
    
    % Build the matrix T(f(x,y))
    T1 = BuildT1(fxy,n1-m1,n2-m2);
     
    arr_T1{i-1} = T1;
end

C = blkdiag(arr_T1{:});

% %
% %
% Form the Right hand side vector

for i = 1:1:nPolys_fxy - 1
    
    % temporarily label the ith entry of the array as fxy
    fxy = GetAsVector(arr_fxy{i});
    
    % Get degree structure of fxy
    [m1,m2] = GetDegree(arr_fxy{i});
        
    arr_rhs{i} = fxy;
end

vRHS = cell2mat(arr_rhs')

% %
% %
% %
x_ls = SolveAx_b(C,vRHS);

% Split solution vector into polynomials h{i}.
nPolys_hxy = nPolys_fxy -1 ;

for i = 1:nPolys_hxy
    
    arr_hxy{i} = x_ls(1:nNonZeros_h{i});
    
    x_ls(1:nNonZeros_h{i}) = [];
    
end
end