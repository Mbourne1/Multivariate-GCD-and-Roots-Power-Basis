
function [fxy_matrix, m] = Examples_Roots_FromCoefficients(ex_num)


x = sym('x');
y = sym('y');


switch ex_num
    case '1'
        f = (x+1) * (y-1) * (y-3);
        
    case '2'
        f = (x + y + 1)^2 * (x+1)^2 * (y-1)^3 * (y-3);
        
    case '3'
        f = (x+0.5)^3 * (y-0.27535)^2 * (x-1.5276543);
        
    case '4'
        f = (x-0.123456)^4 * (x-0.6789)^3 * (x + y - 0.2468)^2 * (y-5.2479)^3;
    otherwise
        error([mfilename 'Not a valid example number'])
end

display(f);
m = double(feval(symengine, 'degree', f));
m1 = double(feval(symengine, 'degree', f, x));
m2 = double(feval(symengine, 'degree', f, y));

fxy_matrix = double(rot90(coeffs(f,[x,y],'All'),2));


end