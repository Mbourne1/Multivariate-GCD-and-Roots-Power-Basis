
function [fxy_matrix, m] = Examples_Roots_FromCoefficients(ex_num)

switch ex_num
    case '1'
        x = sym('x');
        y = sym('y');

        f = (x^2 + y^2 + 1)^2 * (x+1) * (y-1) * (y-3);
        
        m = double(feval(symengine, 'degree', f));
        
        fxy_matrix = double(rot90(coeffs(f,[x,y],'All'),2));
        
    case '2'
        x = sym('x');
        y = sym('y');

        f = (x + y + 1)^2 * (x+1)^2 * (y-1)^3 * (y-3);
        
        m = double(feval(symengine, 'degree', f));
        
        fxy_matrix = double(rot90(coeffs(f,[x,y],'All'),2)); 
end

end