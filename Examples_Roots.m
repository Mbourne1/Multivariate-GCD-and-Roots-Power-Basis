function [fxy_matrix, m] = Examples_Roots(ex)
% Given an example number, get a set of roots, Build the polynomial and
% output the matrix of coefficients.
%
% Inputs.
% 
%
% ex : Example Number
%
%
% Outputs.
%
%
% fxy_matrix : Coefficient matrix of polynomial f(x,y)
%
% m : Total degree of polynomial f(x,y)

polys_x = [];
polys_y = [];
polys_xy = [];

switch ex
    
    case '1'
        % Example 1
        % in File Bivariate Root Finding - By example.tex
        %
        % (x-1)^2
        % (x-3)
        % (y+2)^2
        % (y-1)
        
        roots_f_x =[...
            1  2;
            3  1;
            ];
        
        roots_f_y =[...
            -2  2;
            1   1;
            ];
        
        m = 6;
        

    case '2'
%         ( x - 0.7526)^{2}
%         ( x - 0.1236)^{3}
%         ( x - 1.1000)^{4}
%         ( y - 0.1459)^{2}
%         ( y - 0.9876)^{1}
%         (-0.75 x^{0}y^{0} + 1.00 x^{1}y^{0} +1.00 x^{0}y^{1} +0.00 x^{1}y^{1} )^{1}
%         ( 2.62 x^{0}y^{0} + 1.00 x^{1}y^{0} +1.00 x^{0}y^{1} +0.00 x^{1}y^{1} )^{2}
       
    
        polys_xy{1,1} = [...
            -0.753  1;
            1       0;
        ];
        polys_xy{2,1} = [...
            2.6235  1;
            1       0;
        ];
        polys_xy{3,1} = [...
            2.6235  1;
            1       0;
        ];
        
        roots_f_x = [...
            0.7526  2;
            0.1236  3;
            1.1000  4;
        ];
    
        roots_f_y = [...
            0.1459  2;
            0.9876  1;
        ];
        
        m = 15;
        
end

% Each row of polys x is a simple polynomial in p(x)

fprintf('Roots with respect to x:\n\n')
disp(roots_f_x)
fprintf('Roots with respect to y:\n\n')
disp(roots_f_y)
fprintf('Non Separable roots \n\n')
disp(polys_xy)

polys_x = mult_roots_x(roots_f_x);
polys_y = mult_roots_y(roots_f_y);


polys_f = [polys_x; polys_y; polys_xy];

PrintRoots_Bivariate(roots_f_x,roots_f_y,polys_xy);

fxy_matrix = BuildPoly_NonSeparable(polys_f);


end


