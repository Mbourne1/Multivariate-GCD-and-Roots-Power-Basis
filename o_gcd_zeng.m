function [u,v,w] = o_gcd_zeng(fxy,gxy)
% Prepare the polynomials f(x,y) and g(x,y) for input into Zengs method.
% Note that the inputs must be polynomials in string form

x = sym('x');
y = sym('y');

%f = GetSymbolicFromCoefficients(fxy,x,y);
%g = GetSymbolicFromCoefficients(gxy,x,y);

f = ((x+1)*(x+4)*(x+5));
g = ((x+1)*(x+4)*(x+6));

f = char(expand(f))
g = char(expand(g))


%f = '-45*x*y - 15*x^3*y - 20*x*y^2 + 27*x*y^3 + 9*x^3*y^3 + 12*x*y^4';
%g = '45*x^2*y^2 + 15*x^2*y^3 - 27*x^2*y^4 - 9*x^2*y^5';

[u,v,w] = PolynomialGCD(f,g)

end