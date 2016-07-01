function [fxy_exact, gxy_exact,...
    uxy_exact,vxy_exact,...
    dxy_exact,...
    m,m1,m2,...
    n,n1,n2,...
    t_exact,t1_exact,t2_exact] = Examples_GCD_FromCoefficients(ex_num)

uxy_exact = [];
vxy_exact = [];

syms x y


switch ex_num
    
    case '1'
        f = -45*x*y - 15*x^3*y - 20*x*y^2 + 27*x*y^3 + 9*x^3*y^3 + 12*x*y^4;
        g = 45*x^2*y^2 + 15*x^2*y^3 - 27*x^2*y^4 - 9*x^2*y^5;
        d = -70*x*y + 42*x*y^3;
        
    case '2'
        
        f = (x^2 + y^2 + 1)^2 * (x+1) * (y-1) * (y-3);
        g = (x^2 + y^2 + 1)^2 * (x+1) * (y-2) * (y+1.234)^2;
        d = (x^2 + y^2 + 1)^2 * (x+1);
        
    case '3'
        % Example 1 from "Computing the Greatest common divisor of ...
        % multivariate polynomials over finite fields" Yang
        
        f = (x +  y + 11) * (3*x + y + 1);
        g = (x + y + 4) * (3*x + y + 1);
        d = (3*x + y + 1);
        
        
    case '4'
        
        f = (x + 11)^2 * (y + 1)^3 * (x + 0.6)^3 * (y- 1.6)^1;
        g = (x + 11)^2 * (y + 1)^3 * (x + 1.9)^7 * (y-0.625)^3;
        d = (x + 11)^2 * (y + 1)^3;
        
    case '5'
        
        f = (x^2 + y^2 + 1)^2 * (x+1) * (y-1) * (y-3)^2 * (x+y-0.1534);
        g = (x^2 + y^2 + 1)^2 * (x+1) * (y-2) * (y+1.234)^2;
        d = (x^2 + y^2 + 1)^2 * (x+1); 
        
    case '6'
        f = (x^2 + y^2 + 0.1575)^2 * (x+0.59)^3 * (x+0.8) * (y-1) * (y-3)^2 * (x+y-0.1534);
        g = (x^2 + y^2 + 0.1575)^2 * (x+0.59)^3 * (x + 0.75) * (y-2.75) * (y+1.234)^2;
        d = (x^2 + y^2 + 0.1575)^2 * (x+0.59)^3; 
        
    case '7'
        f = (x^2 + y^2 + 0.1575)^4 * (x + y + 0.575);
        g = (x^2 + y^2 + 0.1575)^4 * (x + y + 0.521);
        d = (x^2 + y^2 + 0.1575)^4 ;
        
    otherwise
        error([mfilename ' : Not a valid example number' ])
end

m = double(feval(symengine, 'degree', f));
n = double(feval(symengine, 'degree', g));
t_exact = double(feval(symengine, 'degree', d));

fxy_exact = double(rot90(coeffs(f,[x,y],'All'),2));
gxy_exact = double(rot90(coeffs(g,[x,y],'All'),2));
dxy_exact = double(rot90(coeffs(d,[x,y],'All'),2));

m1 = double(feval(symengine, 'degree', f,x));
m2 = double(feval(symengine, 'degree', f,y));
n1 = double(feval(symengine, 'degree', g,x));
n2 = double(feval(symengine, 'degree', g,y));
t1_exact = double(feval(symengine, 'degree', d,x));
t2_exact = double(feval(symengine, 'degree', d,y));
end