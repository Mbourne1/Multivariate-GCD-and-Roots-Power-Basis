
function fxy_sym = GetSymBivariatePoly(fxy_matrix)
% Given the coefficients of a bivariate polynomial f(x,y), get the symbolic
% polynomial.

[nRows,nCols] = size(fxy_matrix)

m1 = nRows - 1;
m2 = nCols - 1;

x = sym('x');
y = sym('y');

fxy_sym = 0;

for i = 1:1:nRows
    for j = 1:1:nCols
        fxy_sym = fxy_sym + (fxy_matrix(i,j) * (x^(i-1)) * (y^(j-1)));
    end
end


end