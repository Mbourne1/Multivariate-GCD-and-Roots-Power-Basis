function arr_hxy = Deconvolve_Bivariate_Batch_Respective(arr_fxy)
% Perform the batch deconvolution
%
% Inputs.
%
% arr_fxy : Array of polynomials f(x,y)
%
%
% % Outputs
%
% arr_hxy : (Array of Matrices) Each cell contains a matrix of coefficients
% of polynomial f(x,y)



% Form the left hand side matrix

% Get number of polynomials in the array arr_fxy
nPolys_arr_fxy = size(arr_fxy,1);

% For each of the polynomials excluding the first.
vNCoefficients_hxy = zeros(nPolys_arr_fxy-1,1);

vDeg_x_fxy = zeros(nPolys_arr_fxy, 1);
vDeg_y_fxy = zeros(nPolys_arr_fxy, 1);

for i = 1:1:nPolys_arr_fxy
    [vDeg_x_fxy(i), vDeg_y_fxy(i)] = GetDegree_Bivariate(arr_fxy{i});
end

vDeg_x_hxy = diff(vDeg_x_fxy);
vDeg_y_hxy = diff(vDeg_y_fxy);


v_n1 = zeros(nPolys_arr_fxy - 1, 1);
v_n2 = zeros(nPolys_arr_fxy - 1, 1);

arr_T1 = cell(nPolys_arr_fxy - 1);

for i = 2:1:nPolys_arr_fxy
    
    % Get degree of f{i-1} (The corresponding vector in the RHS)
    [prev_m1, prev_m2] = GetDegree_Bivariate(arr_fxy{i-1});
    
    % Temporarily call the ith entry f(x,y)
    fxy = arr_fxy{i};
    
    % Get degree of f(x,y) with respect to x and y
    [m1, m2] = GetDegree_Bivariate(fxy);
    
    % Get the degree of h{i-1}(x,y)
    v_n1(i-1) = prev_m1 - m1;
    v_n2(i-1) = prev_m2 - m2;

    
    % Get number of coefficients in h_{i}(x,y)
    vNCoefficients_hxy(i-1) = (v_n1(i-1)+1) * (v_n2(i-1)+1);
    
    % Build the matrix T(f(x,y))
    T1 = BuildT1_Relative_Bivariate(fxy, prev_m1-m1, prev_m2-m2);
     
    arr_T1{i-1} = T1;
end

C = blkdiag(arr_T1{:});

% %
% %
% Form the Right hand side vector
vRHS = GetRHS_Vector(arr_fxy);

% %
% %
% %
x_ls = SolveAx_b(C,vRHS);

% Split solution vector into polynomials h{i}.
nPolys_hxy = nPolys_arr_fxy -1 ;

arr_hxy = cell(nPolys_hxy,1);
for i = 1:nPolys_hxy
    
    
    temp_vec = x_ls(1:vNCoefficients_hxy(i));
    
    arr_hxy{i,1} = GetAsMatrix(temp_vec,v_n1(i),v_n2(i));
    
    
    x_ls(1:vNCoefficients_hxy(i)) = [];    
end

end


function vRHS = GetRHS_Vector(arr_fxy)
%
% % Inputs
%
% arr_fxy : (Array of Matrices) Each cell contains matrix of coefficients
% of polynomail f_{i}(x,y)
%
% % Outputs
%
% vRHS : (Vector) Vector containing coefficients of all polynomials
% f_{0}(x,y),...,f_{n-1}(x,y)


% Get number of polynomials in the array 
nPolys_arr_fxy = length(arr_fxy);

for i = 1:1:nPolys_arr_fxy - 1
    
    % temporarily label the ith entry of the array as fxy
    fxy = GetAsVector(arr_fxy{i});
    
    % Get degree structure of fxy
    [m1, m2] = GetDegree_Bivariate(arr_fxy{i});
        
    arr_rhs{i} = fxy;
end

vRHS = cell2mat(arr_rhs');


end