
function fxy_sym = GetSymBivariatePoly(fxy_matrix)
% Given the coefficients of a bivariate polynomial f(x,y), get the symbolic
% polynomial.

% Get number of rows and columns in f(x,y)
[nRows,nCols] = size(fxy_matrix);

% Initialise the symbolic variables x and y
x = sym('x');
y = sym('y');

fxy_sym = 0;

for i = 1:1:nRows
    for j = 1:1:nCols
        fxy_sym = fxy_sym + (fxy_matrix(i,j) * (x^(i-1)) * (y^(j-1)));
    end
end


end