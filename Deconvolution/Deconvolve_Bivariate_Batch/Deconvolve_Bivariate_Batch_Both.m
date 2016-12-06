function arr_hxy = Deconvolve_Bivariate_Batch_Both(arr_fxy, vDeg_t_arr_fxy, vDeg_x_arr_fxy, vDeg_y_arr_fxy)
% Perform the batch deconvolution with the constraints on the total degree
% and relative degree
%
% Inputs.
%
% arr_fxy : Array of polynomials f(x,y)
%
% vDegt_arr_fxy : Vector containing the total degree of the polynomials f_{i}
%
% vDegx_arr_fxy : Vector containing the degree of polynomials f_{i} with
% respect to x
%
% vDegy_arr_fxy : Vector containing the degree of polynomials f_{i} with
% respect to y


% %
% %
% Form the left hand side matrix

% Get the number of polynomials in the array arr_f_{i}(x,y)
nPolys_arr_fxy = size(arr_fxy,1);

% Get the number of polynomials in the array arr_h_{i}(x,y)
nPolys_arr_hxy = nPolys_arr_fxy - 1;

% Get the total degree of polynomials in the array h_{i}(x,y)
vDeg_t_arr_hxy = vDeg_t_arr_fxy(1:end-1) - vDeg_t_arr_fxy(2:end);

% Get the degree of polynomials h_{i}(x,y) with respect to x
vDeg_x_arr_hxy = vDeg_x_arr_fxy(1:end-1) - vDeg_x_arr_fxy(2:end);

% Get the degree of polynomials h_{i}(x,y) with respect to y
vDeg_y_arr_hxy = vDeg_y_arr_fxy(1:end-1) - vDeg_y_arr_fxy(2:end);


% For each of the polynomials excluding the first.
for i = 2:1:nPolys_arr_fxy
    
    % Get total degree of f{i-1} (The corresponding vector in the RHS)
    m_prev = vDeg_t_arr_fxy(i-1);
    
    % Get the degree of f_{i-1} with respect to x
    m1_prev = vDeg_x_arr_fxy(i-1);
    
    % Get the degree of f_{i-1} with respect to y
    m2_prev = vDeg_y_arr_fxy(i-1);
    
    % Temporarily call the ith entry f(x,y)
    fxy = arr_fxy{i};
    
    % Get the total degree
    m_current = vDeg_t_arr_fxy(i);
    
    % Get the degree of f_{i} with respect to x
    m1_current = vDeg_x_arr_fxy(i);
    
    % Get the degree of f_{i} with respect to y
    m2_current = vDeg_y_arr_fxy(i);
    
    n = m_prev - m_current;
    n1 = m1_prev - m1_current;
    n2 = m2_prev - m2_current;
    
    % Build the matrix T(f(x,y))
    T1 = BuildT1_Both(fxy, m_current, n, n1, n2);
    
    % Add T1 to array of matrices
    arr_T1{i-1} = T1;
end

% Build a block diagonal of the matrices
C_fww = blkdiag(arr_T1{:});

% %
% %
% Form the Right hand side vector
arr_rhs = cell(nPolys_arr_fxy - 1,1);

for i = 1:1:nPolys_arr_fxy - 1
    
    % temporarily label the ith entry of the array as fxy
    v_fxy = GetAsVector(arr_fxy{i});
    
    % Get degree structure of fxy
    m = vDeg_t_arr_fxy(i);
    
    % Get degree of f_{i} with respect to x
    m1 = vDeg_x_arr_fxy(i);
    
    % Get degree of f_{i} with respect to y
    m2 = vDeg_y_arr_fxy(i);
    
    % Remove zeros from the end of the vector
    nNonZeros_fxy = GetNumNonZeros(m1,m2,m);
    
    % Remove the zero entries 
    v_fxy = v_fxy(1:nNonZeros_fxy,1);
    
    % Add the vector to the array
    arr_rhs{i} = v_fxy;
end

vRHS = cell2mat(arr_rhs);

% %
% %
% %
x_ls = pinv(C_fww) * vRHS;
%x_ls = SolveAx_b(C_fww,vRHS);

% Split solution vector into polynomials h{i}.
nPolys_arr_hxy = nPolys_arr_fxy - 1 ;

arr_hxy = cell(nPolys_arr_hxy,1);

for i = 1:nPolys_arr_hxy
    
    % Get coefficients as a matrix
    % error('To do')
    
    % Ge the total degree of h_{i}(x,y) 
    n = vDeg_t_arr_hxy(i);
    
    % Get degree of h_{i}(x,y) with respect to x
    n1 = vDeg_x_arr_hxy(i);
    
    % Get degree of h_{i}(x,y) with respect to y
    n2 = vDeg_y_arr_hxy(i);
    
    nNonZeros_hxy = GetNumNonZeros(n1,n2,n);
    
    temp_vec = x_ls(1:nNonZeros_hxy);
    x_ls(1:nNonZeros_hxy) = [];
    
    
    % Append zeros so that temp_vec fills a n1+1,n2+1 matrix
    nCoefficients = (n1+1) * (n2+1);
    nZeros = nCoefficients - nNonZeros_hxy;
    temp_vec = [temp_vec ; zeros(nZeros,1)];
    
    arr_hxy{i} = GetAsMatrix(temp_vec,n1,n2);
    
    
    
    
end
end