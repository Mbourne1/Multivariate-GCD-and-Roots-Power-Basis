function [dist] = GetDistanceBetweenPolynomials(fxy_exact,fxy_computed,name)
% Given the two polynomials f(x) exact and f(x) computed. Calculate the
% distance between them.

LineBreakLarge()

% Ensure that f(x,y) exact and f(x,y) computed are both in term of
% total degree.

fxy_computed = Normalise(fxy_computed);
fxy_exact = Normalise(fxy_exact);

dist = norm(fxy_exact - fxy_computed) ./ norm(fxy_exact);
fprintf([mfilename ' : ' sprintf('Distance of %s exact from %s computed a - b / a : \t %2.4e \n',name,name,dist)]);





%display(fxy_exact./fxy_exact(1,1));
%display(fxy_computed./fxy_computed(1,1));

% Normalise f(x,y) exact
fxy_exact = NormaliseMatrix(fxy_exact);

% Normalise f(x,y) computed
fxy_computed = NormaliseMatrix(fxy_computed);





end

function costheta = GetAngle(a,b)
costheta = dot(a,b)/(norm(a)*norm(b));
end