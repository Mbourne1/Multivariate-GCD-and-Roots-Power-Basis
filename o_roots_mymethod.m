function [] = o_roots_mymethod(fxy_matrix,M)
% Get the factorization of a bivariate polynomial f(x,y)
%
% Inputs.
%
% fxy_matrix : Coefficients of the polynomial f(x,y)
%
% M : Total degree of f(x,y)

% Get degree of polynomial f(x,y)
[m1,m2] = GetDegree(fxy_matrix);

% Get factorization with respect to x
if m1 > 0
    [wx,vDegt_x] = o_roots_mymethod_x(fxy_matrix,M);
else
    fprintf([mfilename ' : ''No roots with respect to x \n'])
end

% Get factorization with respect to y
if m2 > 0
    [wy,vDegt_y] = o_roots_mymethod_y(fxy_matrix,M);
else
    % No roots wrt y
    fprintf([mfilename ' : ''No roots with respect to y \n'])
end

try
    [wx,wy,wxy] = o_roots_mymethod_xy(wx,wy,vDegt_x,vDegt_y);
    
    for i = 1:1:length(wx)
        fprintf([mfilename sprintf(' : Roots of x, of degree : %i', i)])
        wx{i}
    end
    for i = 1:1:length(wy)
        fprintf([mfilename sprintf(' : Roots of y, of degree : %i', i)])
        wy{i}
    end
    for i = 1:1:length(wxy)
        fprintf([mfilename sprintf(' : Roots of x and y, of degree : %i', i)])
        wxy{i}
    end
catch err
    fprintf([mfilename err.message])
    
end

% % 
% % 
% %
% %



end