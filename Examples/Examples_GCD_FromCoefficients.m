function [fxy_exact, gxy_exact,...
    uxy_exact,vxy_exact,...
    dxy_exact,...
    m,m1,m2,...
    n,n1,n2,...
    t_exact,t1_exact,t2_exact] = Examples_GCD_FromCoefficients(ex_num)

uxy_exact = [];
vxy_exact = [];

switch ex_num
    case '1'
        %'-45*x*y - 15*x^3*y - 20*x*y^2 + 27*x*y^3 + 9*x^3*y^3 + 12*x*y^4';
        fxy_exact = ...
            [
            0   0   0   0  0
            0 -45 -20 +27 12
            0   0   0   0  0
            0 -15   0   9  0
            ];
        m = 6;
        m1 = 3;
        m2 = 4;
        % 45*x^2*y^2 + 15*x^2*y^3 - 27*x^2*y^4 - 9*x^2*y^5
        gxy_exact = ...
            [
            0 0  0  0   0  0
            0 0  0  0   0  0
            0 0 45 15 -27 -9
            ];
        n = 7;
        n1 = 2;
        n2 = 5;
        
        %'-70*x*y + 42*x*y^3'
        dxy_exact = ...
            [
            0 0 0 0;
            0 -70 0 42
            ]
        t_exact = 4;
        t1_exact = 1;
        t2_exact = 3;
        
    case '2'
        
        fxy_exact = ...
            [
                2 -1
                -3 0
                1 0
            ];
        m = 2;
        m1 = 2;
        m2 = 1;
        
        gxy_exact = ...
            [
                3   -1
                -4  0
                1   0
            ];
        n = 2;
        n1 = 2;
        n2 = 1;
        
        dxy_exact = ...
            [
            -1  
            1
            ];
        t_exact = 1;
        t1_exact = 1;
        t2_exact = 0;
        
    case '3'
        x = sym('x');
        y = sym('y');

        f = (x^2 + y^2 + 1)^2 * (x+1) * (y-1) * (y-3);
        g = (x^2 + y^2 + 1)^2 * (x+1) * (y-2);
        d = (x^2 + y^2 + 1)^2 * (x+1);
        
        m = double(feval(symengine, 'degree', f));
        n = double(feval(symengine, 'degree', g));
        t_exact = double(feval(symengine, 'degree', d));
        
        fxy_exact = double(rot90(coeffs(f,[x,y],'All'),2));
        gxy_exact = double(rot90(coeffs(g,[x,y],'All'),2));
        dxy_exact = double(rot90(coeffs(d,[x,y],'All'),2));
        
        m1 = 5;
        m2 = 6;
        n1 = 5;
        n2 = 5;
        t1_exact = 5;
        t2_exact = 4;
        
    case '4'
        % Example 1 from "Computing the Greatest common divisor of ...
        % multivariate polynomials over finite fields" Yang 
        
        x = sym('x');
        y = sym('y');

        f = (x +  y + 11) * (3*x + y + 1);
        g = (x + y + 4) * (3*x + y + 1);
        d = (3*x + y + 1);
        
        m = double(feval(symengine, 'degree', f));
        n = double(feval(symengine, 'degree', g));
        t_exact = double(feval(symengine, 'degree', d));
        
        fxy_exact = double(rot90(coeffs(f,[x,y],'All'),2));
        gxy_exact = double(rot90(coeffs(g,[x,y],'All'),2));
        dxy_exact = double(rot90(coeffs(d,[x,y],'All'),2));
        
        m1= 0;
        m2 = 0;
        n1 = 0;
        n2 = 0;
        t1_exact = 0;
        t2_exact = 0;
    case '5'
        
        x = sym('x');
        y = sym('y');

        f = (x + 11)^2 * (y + 1)^3 * (x + 0.6)^3 * (y- 1.6)^1;
        g = (x + 11)^2 * (y + 1)^3 * (x + 1.9)^7 * (y-0.625)^3;
        d = (x + 11)^2 * (y + 1)^3;
        
        m = double(feval(symengine, 'degree', f));
        n = double(feval(symengine, 'degree', g));
        t_exact = double(feval(symengine, 'degree', d));
        
        fxy_exact = double(rot90(coeffs(f,[x,y],'All'),2));
        gxy_exact = double(rot90(coeffs(g,[x,y],'All'),2));
        dxy_exact = double(rot90(coeffs(d,[x,y],'All'),2));
        
        m1= 0;
        m2 = 0;
        n1 = 0;
        n2 = 0;
        t1_exact = 0;
        t2_exact = 0;
end
end