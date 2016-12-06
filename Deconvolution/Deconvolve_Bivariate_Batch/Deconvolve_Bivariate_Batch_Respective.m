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
nPolys_arr_fxy = size(arr_fxy,1);

% For each of the polynomials excluding the first.
vNumCoefficients_hxy = zeros(nPolys_arr_fxy-1,1);
v_n1 = zeros(nPolys_arr_fxy-1,1);
v_n2 = zeros(nPolys_arr_fxy-1,1);

for i = 2:1:nPolys_arr_fxy
    
    % Get degree of f{i-1} (The corresponding vector in the RHS)
    [prev_m1,prev_m2] = GetDegree(arr_fxy{i-1});
    
    % Temporarily call the ith entry f(x,y)
    fxy = arr_fxy{i};
    
    % Get degree of f(x,y) with respect to x and y
    [m1,m2] = GetDegree(fxy);
    
    % Get the degree of h{i-1}(x,y)
    v_n1(i-1) = prev_m1 - m1;
    v_n2(i-1) = prev_m2 - m2;

    
    % Get number of coefficients in h_{i}(x,y)
    vNumCoefficients_hxy(i-1) = (v_n1(i-1)+1) * (v_n2(i-1)+1);
    
    % Build the matrix T(f(x,y))
    T1 = BuildT1_Relative(fxy,prev_m1-m1,prev_m2-m2);
     
    arr_T1{i-1} = T1;
end

C = blkdiag(arr_T1{:});

% %
% %
% Form the Right hand side vector

for i = 1:1:nPolys_arr_fxy - 1
    
    % temporarily label the ith entry of the array as fxy
    fxy = GetAsVector(arr_fxy{i});
    
    % Get degree structure of fxy
    [m1,m2] = GetDegree(arr_fxy{i});
        
    arr_rhs{i} = fxy;
end

vRHS = cell2mat(arr_rhs');

% %
% %
% %
x_ls = SolveAx_b(C,vRHS);

% Split solution vector into polynomials h{i}.
nPolys_hxy = nPolys_arr_fxy -1 ;

for i = 1:nPolys_hxy
    
    
    temp_vec = x_ls(1:vNumCoefficients_hxy(i));
    
    arr_hxy{i,1} = GetAsMatrix(temp_vec,v_n1(i),v_n2(i));
    
    
    x_ls(1:vNumCoefficients_hxy(i)) = [];    
end
end