function [fxy_matrix] = BuildPoly_NonSeparable(roots)
% Given a set of bivariate polynomial roots where each root is in matrix
% form, obtain the coefficients of the original polynomial.

% Get the number of roots in the set
num_roots = size(roots,1);

% initalise process by setting the first root to be the roots to be
% multiplied by.
new_p1_mtrx = roots{1,1};

% for all remaining roots, in turn, multiply
for i = 2:1:num_roots
    
    % let p2 be the polynomial of the root which is being multiplied by p1
    p2 = roots{i,1};
    
    % get dimensions of the poly so far
    [r,c] = size(new_p1_mtrx);
    
    
    % Get degree of the convolved polynomial with respect to x and y
    deg_p1_x = r-1;
    deg_p1_y = c-1;
    
    % Get degree of the new polynomial to be convolved
    [p2_r, p2_c] = size(p2);
    deg_p2_x = p2_r - 1;
    deg_p2_y = p2_c - 1;
    
    % % form a vector from the coefficients of polynomial p1
    % tot = num of diagonals in poly
    vec_poly_p1 = [];
    tot = r + c - 1;
    for t = 0:1:tot
        for i2 = 0:1:t
            i1 = t -i2;
            if (i1< r) && (i2 < c)
                vec_poly_p1 = [vec_poly_p1 ; new_p1_mtrx(i1+1,i2+1)];
            end
        end
    end
    
    % Get the new degrees of polynomial p1, when multiplied by p2
    new_deg_p1_x = deg_p1_x + deg_p2_x;
    new_deg_p1_y = deg_p1_y + deg_p2_y;
    
    T1 = BuildT1(p2,r-1,c-1);
    
    
    new_p1_vec = T1 * (vec_poly_p1);
    
    new_p1_vec = new_p1_vec';
    
    % rebuild polynomial 1 as a matrix
    rows = new_deg_p1_x + 1;
    cols = new_deg_p1_y + 1;
    
    diags = rows + cols - 1;
    
    count = 1;
    for tot = 0:1:diags-1
        for i3 = tot:-1:0
            j3 = tot - i3;
            if i3<rows && j3<cols
                new_p1_mtrx(i3+1,j3+1) =  new_p1_vec(count);
                count = count + 1;
            end
        end
    end
    
end

fxy_matrix = new_p1_mtrx;

end