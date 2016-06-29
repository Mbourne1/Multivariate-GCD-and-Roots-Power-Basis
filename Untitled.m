function [] = Untitled()

% f = x^5 y-x^5+x^4 y-x^4+2 x^3 y^3-2 x^3 y^2+2 x^3 y-2 x^3+2 x^2 y^3-2 x^2 y^2+2 x^2 y-2 x^2+x y^5-x y^4+2 x y^3-2 x y^2+x y-x+y^5-y^4+2 y^3-2 y^2+y-1
% f = (x^2 + y^2 + 1)^2 (x+1)(y-1)
% %
% g = (x^2 + y^2 + 1)^2 (x+1) (y-2)
% g = x^5 y-2 x^5+x^4 y-2 x^4+2 x^3 y^3-4 x^3 y^2+2 x^3 y-4 x^3+2 x^2 y^3-4 x^2 y^2+2 x^2 y-4 x^2+x y^5-2 x y^4+2 x y^3-4 x y^2+x y-2 x+y^5-2 y^4+2 y^3-4 y^2+y-2

x = sym('x');
y = sym('y');

f = (x^2 + y^2 + 1)^2 * (x+1) * (y-1) * (y-3);
g = (x^2 + y^2 + 1)^2 * (x+1) * (y-2);

m = feval(symengine, 'degree', f);
n = feval(symengine, 'degree', g);

f_mat = rot90(coeffs(f,[x,y],'All'),2);
g_mat = rot90(coeffs(g,[x,y],'All'),2);


limits_t = [0,min(m,n)];

o_gcd_mymethod(f_mat,g_mat,m,n,limits_t);

end