function [fxy, gxy,uxy,vxy,dxy,...
    m,m1,m2,...
    n,n1,n2,...
    t,t1,t2] = Examples_GCD_FromCoefficients(ex_num)

uxy = [];
vxy = [];

syms x y

addpath('../Examples')
[f,g,d,u,v] = Bivariate_GCD_Examples(ex_num);

symbolic_d = GetSymbolicPoly(d);
symbolic_f = GetSymbolicPoly(f);
symbolic_g = GetSymbolicPoly(g);

% Get coefficients 
fxy = double(rot90(coeffs(symbolic_f,[x,y],'All'),2));
gxy = double(rot90(coeffs(symbolic_g,[x,y],'All'),2));
dxy = double(rot90(coeffs(symbolic_d,[x,y],'All'),2));


% Get degree of f(x,y), g(x,y) and d(x,y)
m = double(feval(symengine, 'degree', symbolic_f));
n = double(feval(symengine, 'degree', symbolic_g));
t = double(feval(symengine, 'degree', symbolic_d));

display(symbolic_f);
display(symbolic_g);
display(symbolic_d);

m = double(feval(symengine, 'degree', symbolic_f));
n = double(feval(symengine, 'degree', symbolic_g));
t = double(feval(symengine, 'degree', symbolic_d));



m1 = double(feval(symengine, 'degree', symbolic_f,x));
m2 = double(feval(symengine, 'degree', symbolic_f,y));
n1 = double(feval(symengine, 'degree', symbolic_g,x));
n2 = double(feval(symengine, 'degree', symbolic_g,y));
t1 = double(feval(symengine, 'degree', symbolic_d,x));
t2 = double(feval(symengine, 'degree', symbolic_d,y));



end