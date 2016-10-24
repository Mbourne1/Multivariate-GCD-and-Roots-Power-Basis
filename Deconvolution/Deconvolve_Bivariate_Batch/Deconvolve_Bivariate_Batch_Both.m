function arr_hxy = Deconvolve_Bivariate_Batch_Both(arr_fxy,vDegt_arr_fxy)
% Perform the batch deconvolution with the constraints on the total degree
% and relative degree
%
% Inputs.
%
% arr_fxy : Array of polynomials f(x,y)
%
% vDegt_arr_fxy : Vector containing the total degree of the polynomials f{i}

% %
% %
% Form the left hand side matrix

% Get number of polynomials in the array arr_fxy
nPolys_arr_fxy = size(arr_fxy,1);
nPolys_arr_hxy = nPolys_arr_fxy - 1;

nNonZeros_hxy = zeros(nPolys_arr_hxy,1);

vDegt1_arr_fxy = zeros(nPolys_arr_fxy,1);
vDegt2_arr_fxy = zeros(nPolys_arr_fxy,1);

for i = 1:1:nPolys_arr_fxy
    
    % Get degree of polynomial f_{i}(x,y) with respect to x and to y.
    [vDegt1_arr_fxy(i,1),vDegt2_arr_fxy(i,1)] = GetDegree(arr_fxy{i});
   
end

vDegt_arr_hxy = vDegt_arr_fxy(1:end-1) - vDegt_arr_fxy(2:end);
vDegt1_arr_hxy = vDegt1_arr_fxy(1:end-1) - vDegt1_arr_fxy(2:end);
vDegt2_arr_hxy = vDegt2_arr_fxy(1:end-1) - vDegt2_arr_fxy(2:end);


% For each of the polynomials excluding the first.
for i = 2:1:nPolys_arr_fxy
    
    % Get degree of f{i-1} (The corresponding vector in the RHS)
    m_prev = vDegt_arr_fxy(i-1);
    [m1_prev,m2_prev] = GetDegree(arr_fxy{i-1});
    
    % Temporarily call the ith entry f(x,y)
    fxy = arr_fxy{i};
    
    % Build the matrix T(f(x,y))
    T1 = BuildT1(fxy,vDegt1_arr_hxy(i-1),vDegt2_arr_hxy(i-1));
    
    % % Remove the columns corresponding to zeros in f{i-1}
    
    % Get number of zeros in h{i-1}
    nNonZeros_hxy(i-1) = GetNumNonZeros(vDegt1_arr_hxy(i-1),vDegt2_arr_hxy(i-1),vDegt_arr_hxy(i-1));
    
    % Remove the columns
    T1 = T1(:,1:nNonZeros_hxy(i-1));
    
    % % Remove the rows which correspond to zeros in f{i-1}
    nNonZeros = GetNumNonZeros(m1_prev,m2_prev,m_prev);
    
    % Remove the rows
    T1 = T1(1:nNonZeros,:);
    
    arr_T1{i-1} = T1;
end

C_fww = blkdiag(arr_T1{:});

% %
% %
% Form the Right hand side vector
arr_rhs = cell(nPolys_arr_fxy - 1);

for i = 1:1:nPolys_arr_fxy - 1
    
    % temporarily label the ith entry of the array as fxy
    fxy = GetAsVector(arr_fxy{i});
    
    % Get degree structure of fxy
    m = vDegt_arr_fxy(i);
    [m1,m2] = GetDegree(arr_fxy{i});
    
    % Remove zeros from the end of the vector
    nNonZeros = GetNumNonZeros(m1,m2,m);
    
    % Remove the zero entries
    fxy = fxy(1:nNonZeros,1);
    
    arr_rhs{i} = fxy;
end

vRHS = cell2mat(arr_rhs);

% %
% %
% %
x_ls = SolveAx_b(C_fww,vRHS);

% Split solution vector into polynomials h{i}.
nPolys_hxy = nPolys_arr_fxy -1 ;

arr_hxy = cell(nPolys_hxy,1);
for i = 1:nPolys_hxy
    
    % Get coefficients as a matrix
    % error('To do')
    
    n1 = vDegt1_arr_hxy(i);
    n2 = vDegt2_arr_hxy(i);
    
    temp_vec = x_ls(1:nNonZeros_hxy(i));
    
    % Append zeros so that temp_vec fills a n1+1,n2+1 matrix
    nCoefficients = (n1+1) * (n2+1);
    nZeros = nCoefficients - nNonZeros_hxy(i);
    temp_vec = [temp_vec ; zeros(nZeros,1)];
    
    arr_hxy{i} = GetAsMatrix(temp_vec,n1,n2);
    
    x_ls(1:nNonZeros_hxy(i)) = [];
    
    
end
end