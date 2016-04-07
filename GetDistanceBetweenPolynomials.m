function [dist] = GetDistanceBetweenPolynomials(fxy_exact,fxy_computed,name)
% Given the two polynomials f(x) exact and f(x) computed. Calculate the
% distance between them.

% Normalise f(x,y) exact.
fxy_exact = NormaliseMatrix(fxy_exact);

% Normalise f(x,y) computed.
fxy_computed = NormaliseMatrix(fxy_computed);

dist = norm(fxy_exact - fxy_computed) ./ norm(fxy_exact);

fprintf('Distance of %s exact from %s computed a - b / a : \t %2.4e \n',name,name,dist);
end