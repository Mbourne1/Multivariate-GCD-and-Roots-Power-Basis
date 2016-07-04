function arr_hxy = Deconvolve_Bivariate_Batch_Both(arr_fxy,vDeg_fxy)
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
    n = vDeg_fxy(i-1);
    [n1,n2] = GetDegree(arr_fxy{i-1});
    
    % Temporarily call the ith entry f(x,y)
    fxy = arr_fxy{i};
    
    % Get the total degree of f(x,y) from the input vector
    m = vDeg_fxy(i);
    
    % Get degree of f(x,y) with respect to x and y
    [m1,m2] = GetDegree(fxy);
    
    % Build the matrix T(f(x,y))
    T1 = BuildT1(fxy,n1-m1,n2-m2);
    
    % % Strip the columns corresponding to zeros in f{i-1}
    
    % Get number of zeros in h{i-1}
    nNonZeros_h{i-1} = GetNumNonZeros(n1-m1,n2-m2,n-m);
    
    % Remove the columns
    T1 = T1(:,1:nNonZeros_h{i-1});
    
    % Remove the rows which correspond to zeros in f{i-1}
    nNonZeros = GetNumNonZeros(n1,n2,n);
    
    % Remove the rows
    T1 = T1(1:nNonZeros,:);
    
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
    m = vDeg_fxy(i);
    [m1,m2] = GetDegree(arr_fxy{i});
    
    % Remove zeros from the end of the vector
    nNonZeros = GetNumNonZeros(m1,m2,m);
    
    % Remove the zero entries
    fxy = fxy(1:nNonZeros,1);
    
    arr_rhs{i} = fxy;
end

vRHS = cell2mat(arr_rhs');

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