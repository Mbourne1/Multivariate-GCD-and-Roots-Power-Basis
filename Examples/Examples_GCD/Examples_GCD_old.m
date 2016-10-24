function [fxy_mtrx_exct,gxy_mtrx_exct,...
    uxy_mtrx_exct, vxy_mtrx_exct,...
    dxy_mtrx_exct,...
    m,m1,m2,...
    n,n1,n2,...
    t,t1,t2] = Examples_GCD(ex)

% NOTE -THESE EXAMPLES ASSUME SMALLEST POWER FIRST


switch ex
       
    case '0'
        % Easy Separable Example
        % Section 1.1
        % f_{x,1} and f'_{x,2}
        
        
        
        fxy_mtrx_exct = ...
            [
            12   0   -9   -3;
            -28   0   21    7;
            20   0  -15   -5;
            -4   0    3    1;
            ];
        
        gxy_mtrx_exct = ...
            [
            -28     0    21       7 ;
            40     0   -30     -10 ;
            -12     0     9       3 ;
            ];
        dxy_mtrx_exct = ...
            [
            4   0   -3  -1  ;
            -4   0    3   1  ;
            ];
        uxy_mtrx_exct = ...
            [
            3   ;
            -4  ;
            1   ;
            ];
        
        vxy_mtrx_exct = ...
            [
            -7;
            3;
            ];
        
        m = 6;
        n = 5;
        t = 4;
        
    case '1'
        % Easy Separable Example
        % f_{x,2} and f'_{x,2}
        % INCORRECT DEGREE CALCULATED
        fxy_mtrx_exct = ...
            [
            1   3   0   -4;
            -1  -3  0   4;
            ];
        gxy_mtrx_exct = ...
            [
            1   3   0   -4;
            ];
        uxy_mtrx_exct = ...
            [
            1   ;
            -1  ;
            ];
        vxy_mtrx_exct = ...
            1;
        
        dxy_mtrx_exct = ...
            [
            1   3   0   -4;
            ];
        m = 4;
        n = 3;
        t = 3;
        
    case '2'
        % Section 1.2
        % Easy Separable Example
        % f_{y,1}(x,y) and f'_{y,1}(x,y) to obtain f_{y,2}
        fxy_mtrx_exct = ...
            [
            12   0   -9   -3;
            -28   0   21    7;
            20   0  -15   -5;
            -4   0    3    1;
            ];
        
        gxy_mtrx_exct = ...
            [
            0   -18     -9;
            0   42      21;
            0   -30     -15;
            0   6       3;
            ];
        
        uxy_mtrx_exct = ...
            [
            2   3   1  ;
            ];
        
        vxy_mtrx_exct = ...
            [
            0   3;
            ];
        
        dxy_mtrx_exct = ...
            [
            -6      -3;
            14       7;
            -10     -5;
            2        1;
            ];
        
        m = 6;
        n = 5;
        t = 4;
        
    case '3'
        % Bivariate - NonSeparable - Easier Example
        % Section 2.1
        %f_{x,1} and f'_{x,1}
        fxy_mtrx_exct = ...
            [
            -5      6      -1   ;
            10     -7       1   ;
            -5      1       0   ;
            ];
        gxy_mtrx_exct = ...
            [
            10  -7  1;
            -10 2   0;
            ];
        uxy_mtrx_exct = ...
            [
            1   -1;
            -2  1;
            1   0;
            ];
        
        vxy_mtrx_exct = ...
            [
            -2  1;
            2   0;
            ];
        
        dxy_mtrx_exct = ...
            [
            -5  1;
            ];
        m = 3;
        n = 2;
        t = 1;
        
        
    case '4'
        % Bivariate Case
        % Section 3.1
        % Non Separable - Difficult Example
        % f_{x,1} and f'_{x,1}
        fprintf('This example relates to Section 3.1 of the Internal Report\n')
        
        fxy_mtrx_exct = ...
            [
            20  -24 -1  6   -1;
            -40 28  6   -7  1;
            15  2   -6  1   0;
            10  -7  1   0   0;
            -5  1   0   0   0;
            ];
        gxy_mtrx_exct = ...
            [
            -40     28      6       -7  1;
            30      4       -12     2   0;
            30      -21     3       0   0;
            -20     4       0       0   0;
            ];
        uxy_mtrx_exct = ...
            [
            4   4   1   1;
            8   4   2   1;
            3   1   1   0;
            -2  1   0   0;
            1   0   0   0;
            ];
        vxy_mtrx_exct = ...
            [
            8   -4  -2  1;
            -6  -2   2  0;
            -6  3   0   0;
            4   0   0   0;
            ];
        dxy_mtrx_exct = ...
            [
            -5  1;
            ];
        m = 5;
        n = 4;
        t = 1;
        
    case '5'
        % Section 3.2 - Partial wrt y
        % Non-Separable - Difficult example
        % GCD of f_{y,1} and f'_{y,1}
        
        % f = (x^2 + y^2 - 4) (x+y-1)(x-1)(y-5)
        % g = f'
        
        fxy_mtrx_exct = ...
            [
            20  -24 -1  6   -1;
            -40 28  6   -7  1;
            15  2   -6  1   0;
            10  -7  1   0   0;
            -5  1   0   0   0;
            ];
        
        gxy_mtrx_exct = ...
            [
            -24   -2      18      -4;
            28    12      -21     4;
            2     -12     3       0;
            -7     2      0       0;
            1     0       0       0;
            ];
        uxy_mtrx_exct = ...
            [
            -20 24  1   -6  1;
            20  -4  -5  1   0;
            5   -6  1   0   0;
            -5  1   0   0   0;
            ];
        vxy_mtrx_exct = ...
            [
            24  2   -18 4;
            -4  -10 3   0;
            -6  2   0   0;
            1   0   0   0
            ];
        dxy_mtrx_exct = ...
            [
            -1 ;
            1 ;
            ];
        
        m = 5 ;
        n = 4 ;
        t = 1 ;
        
        
        
        
        
    case '7'
        % from file anotherexample section 1
        fxy_mtrx_exct = ...
            [
            -2      9       -16     14      -6      1   ;
            7       -29     47      -37     14      -2  ;
            -9      34      -49     33      -10     1   ;
            5       -17     21      -11     2       0   ;
            -1      3       -3      1       0       0   ;
            ];
        gxy_mtrx_exct = ...
            [
            -1  3   -3  1
            ];
        dxy_mtrx_exct = ...
            [
            -1  3   -3  1
            ];
        uxy_mtrx_exct = ...
            [
            2   -3  1;
            -7  8   -2;
            9   -7  1;
            -5  2   0;
            1   0   0;
            ];
        vxy_mtrx_exct = ...
            1;
        
        m = 9;
        n = 3;
        t = 3;
        
    case '8'
        % from file anotherexample section 2
        fxy_mtrx_exct = ...
            [
            -2      9       -16     14      -6      1   ;
            7       -29     47      -37     14      -2  ;
            -9      34      -49     33      -10     1   ;
            5       -17     21      -11     2       0   ;
            -1      3       -3      1       0       0   ;
            ];
        gxy_mtrx_exct = ...
            [
            4    -4  1;
            -2   1   0;
            ];
        dxy_mtrx_exct = ...
            [
            -2  1;
            1   0;
            ];
        uxy_mtrx_exct = ...
            [
            1   -4     6    -4       1;
            -3   11   -15     9      -2;
            3  -10    12     6       1;
            -1    3    -3     1       0;
            ];
        vxy_mtrx_exct = ...
            [
            -2  1;
            ];
        %         m = 7;
        %         n = 2;
        %         t = 1;
        m = 9;
        n = 3;
        t = 1;
        
    case '9'
        % from file anotherexample section 3
        fxy_mtrx_exct = ...
            [
            -2      9       -16     14      -6      1   ;
            7       -29     47      -37     14      -2  ;
            -9      34      -49     33      -10     1   ;
            5       -17     21      -11     2       0   ;
            -1      3       -3      1       0       0   ;
            ];
        gxy_mtrx_exct = ...
            [
            4    -4      1   ;
            -10   9      -2  ;
            8    -6      1   ;
            -2   1       0   ;
            ];
        dxy_mtrx_exct = ...
            [
            -2  1   ;
            5   -2  ;
            -4  1   ;
            1   0   ;
            ];
        uxy_mtrx_exct = ...
            [
            1  -4  6   -4  1;
            -1 3   -3  1   0;
            ];
        vxy_mtrx_exct = ...
            [
            -2  1;
            ];
        
        % 1 - Maximum degree of coefficients
        % 2 - m_{1} + m_{2}
        degree_method = 2;
        switch degree_method
            case 1
                m = 7;
                n = 4;
                t = 3;
            case 2
                m = 9;
                n = 5;
                t = 4;
                
        end
    case '10'
        
        % from file anotherexample section 4
        fxy_mtrx_exct = ...
            [
            -36 -27 13  3   -1;
            9   53  -4  -7  1;
            40  -17 -15 4   0;
            -10 -13 6   0   0;
            -4  4   0   0   0;
            1   0   0   0   0;
            ];
        gxy_mtrx_exct = ...
            [
            -96 112 6   -23 0   1;
            16  28  -47 -1  4   0;
            22  -25 -3  6   0   0;
            -1  -3  4   0   0   0;
            -1  1   0   0   0   0;
            ];
        dxy_mtrx_exct = ...
            [
            -12 -1  1;
            -1  2   0;
            1   0   0;
            ];
        uxy_mtrx_exct = ...
            [
            3  2   -1;
            -1 -4  1;
            -3 2   0;
            1  0   0;
            ];
        vxy_mtrx_exct = ...
            [
            8    -10 1   1;
            -2   0   2   0;
            -1   1   0   0;
            ];
        
        degree_method = 2;
        switch degree_method
            case 1
                m = 5;
                n = 5;
                t = 2;
            case 2
                m = 9;
                n = 9;
                t = 4;
        end
        
    case '11'
        
        % Get roots of f
        roots_f{1,1} = ...
            [
            -1  1;
            1   0;
            ];
        roots_f{2,1} = ...
            [
            -2  1;
            1   0;
            ];
        roots_f{3,1} = ...
            [
            -3  1;
            1   0;
            ];
        fxy_mtrx_exct = BuildPoly_NonSeparable(roots_f);
        
        % Get roots of g
        roots_g{1,1} = ...
            [
            -1  1;
            1   0;
            ];
        roots_g{2,1} = ...
            [
            -2  1;
            1   0;
            ];
        roots_g{3,1} = ...
            [
            -2;
            1;
            ];
        gxy_mtrx_exct = BuildPoly_NonSeparable(roots_g);
        
        % Get roots of polynomial d
        roots_d{1,1} = ...
            [
            -1  1;
            1   0;
            ];
        roots_d{2,1} = ...
            [
            -2  1;
            1   0;
            ];
        
        dxy_mtrx_exct = BuildPoly_NonSeparable(roots_d);
        
        roots_u{1,1} = ...
            [
            -3   1;
            1   0;
            ];
        
        uxy_mtrx_exct = BuildPoly_NonSeparable(roots_u);
        
        roots_v{1,1} = ...
            [
            -2;
            1;
            ];
        
        vxy_mtrx_exct = BuildPoly_NonSeparable(roots_v);
        degree_method = 1;
        switch degree_method
            case 1
                m = 3;
                n = 3;
                t = 2;
            case 2
                m = 6;
                n = 5;
                t = 4;
        end
        
    case '12'
        
        roots_f{1,1} = ...
            [
            -2  0   1;
            0   0   0;
            1   0   0;
            ];
        roots_f{2,1} = ...
            [
            2  0   1;
            0   0   0;
            1   0   0;
            ];
        roots_f{3,1} = ...
            [
            -1;
            1
            ];
        
        roots_g{1,1} = ...
            [
            -2  0   1;
            0   0   0;
            1   0   0;
            ];
        roots_g{2,1} = ...
            [
            2  0   1;
            0   0   0;
            1   0   0;
            ];
        roots_g{3,1} = ...
            [
            -2  1;
            ];
        
        roots_d{1,1} = ...
            [
            -2  0   1;
            0   0   0;
            1   0   0;
            ];
        roots_d{2,1} = ...
            [
            2  0   1;
            0   0   0;
            1   0   0;
            ];
        
        roots_v{1,1} = ...
            [
            -2  1
            ];
        
        roots_u{1,1} = ...
            [
            -1;
            1
            ];
        
        fxy_mtrx_exct = BuildPoly_NonSeparable(roots_f);
        gxy_mtrx_exct = BuildPoly_NonSeparable(roots_g);
        uxy_mtrx_exct = BuildPoly_NonSeparable(roots_u);
        vxy_mtrx_exct = BuildPoly_NonSeparable(roots_v);
        dxy_mtrx_exct = BuildPoly_NonSeparable(roots_d);
        
        m = 5;
        n = 5;
        t = 4;
        
    case '13'
        
        roots_f{1,1} = ...
            [
            -2  0   1;
            0   0   0;
            1   0   0;
            ];
        roots_f{2,1} = ...
            [
            2  0   1;
            0   0   0;
            1   0   0;
            ];
        roots_f{3,1} = ...
            [
            -1;
            1
            ];
        roots_f{4,1} = ...
            [
            -1;
            1
            ];
        
        
        roots_g{1,1} = ...
            [
            -2  0   1;
            0   0   0;
            1   0   0;
            ];
        roots_g{2,1} = ...
            [
            2  0   1;
            0   0   0;
            1   0   0;
            ];
        roots_g{3,1} = ...
            [
            -2  1;
            ];
        roots_g{4,1} = ...
            [
            -5;
            1
            ];
        
        
        roots_d{1,1} = ...
            [
            -2  0   1;
            0   0   0;
            1   0   0;
            ];
        roots_d{2,1} = ...
            [
            2  0   1;
            0   0   0;
            1   0   0;
            ];
        
        roots_v{1,1} = ...
            [
            -2  1
            ];
        roots_v{2,1} = ...
            [
            -5;
            1
            ];
        
        roots_u{1,1} = ...
            [
            -1;
            1
            ];
        roots_u{2,1} = ...
            [
            -1;
            1
            ];
        
        fxy_mtrx_exct = BuildPoly_NonSeparable(roots_f);
        gxy_mtrx_exct = BuildPoly_NonSeparable(roots_g);
        uxy_mtrx_exct = BuildPoly_NonSeparable(roots_u);
        vxy_mtrx_exct = BuildPoly_NonSeparable(roots_v);
        dxy_mtrx_exct = BuildPoly_NonSeparable(roots_d);
        
        m = 6;
        n = 6;
        t = 4;
    case '14'
        % same as example 13 but divided by 10
        
        roots_f{1,1} = ...
            [
            -0.2  0   1;
            0   0   0;
            1   0   0;
            ];
        roots_f{2,1} = ...
            [
            0.2  0   1;
            0   0   0;
            1   0   0;
            ];
        roots_f{3,1} = ...
            [
            -0.1;
            1
            ];
        roots_f{4,1} = ...
            [
            -0.1;
            1
            ];
        
        
        roots_g{1,1} = ...
            [
            -0.2  0   1;
            0   0   0;
            1   0   0;
            ];
        roots_g{2,1} = ...
            [
            0.2  0   1;
            0   0   0;
            1   0   0;
            ];
        roots_g{3,1} = ...
            [
            -0.2  1;
            ];
        roots_g{4,1} = ...
            [
            -0.5;
            1
            ];
        
        
        roots_d{1,1} = ...
            [
            -0.2  0   1;
            0   0   0;
            1   0   0;
            ];
        roots_d{2,1} = ...
            [
            0.2  0   1;
            0   0   0;
            1   0   0;
            ];
        
        roots_v{1,1} = ...
            [
            -0.2  1
            ];
        roots_v{2,1} = ...
            [
            -0.5;
            1
            ];
        
        roots_u{1,1} = ...
            [
            -0.1;
            1
            ];
        roots_u{2,1} = ...
            [
            -0.1;
            1
            ];
        
        fxy_mtrx_exct = BuildPoly_NonSeparable(roots_f);
        gxy_mtrx_exct = BuildPoly_NonSeparable(roots_g);
        uxy_mtrx_exct = BuildPoly_NonSeparable(roots_u);
        vxy_mtrx_exct = BuildPoly_NonSeparable(roots_v);
        dxy_mtrx_exct = BuildPoly_NonSeparable(roots_d);
        
        m = 6;
        n = 6;
        t = 4;
        
        
    case '15'
        % same as example 13 but divided by 10
        
        roots_f{1,1} = ...
            [
            -0.2105  0   1;
            0   0   0;
            1   0   0;
            ];
        roots_f{2,1} = ...
            [
            0.2105  0   1;
            0   0   0;
            1   0   0;
            ];
        roots_f{3,1} = ...
            [
            -0.1;
            1
            ];
        roots_f{4,1} = ...
            [
            -0.1;
            1
            ];
        
        
        roots_g{1,1} = ...
            [
            -0.2105  0   1;
            0   0   0;
            1   0   0;
            ];
        roots_g{2,1} = ...
            [
            0.2105  0   1;
            0   0   0;
            1   0   0;
            ];
        roots_g{3,1} = ...
            [
            -0.2  1;
            ];
        roots_g{4,1} = ...
            [
            -0.5;
            1
            ];
        roots_g{5,1} = ...
            [
            -0.51234;
            1
            ];
        
        
        roots_d{1,1} = ...
            [
            -0.2105  0   1;
            0   0   0;
            1   0   0;
            ];
        roots_d{2,1} = ...
            [
            0.2105  0   1;
            0   0   0;
            1   0   0;
            ];
        
        roots_v{1,1} = ...
            [
            -0.2  1
            ];
        roots_v{2,1} = ...
            [
            -0.5;
            1
            ];
        roots_v{3,1} = ...
            [
            -0.51234;
            1
            ];
        
        roots_u{1,1} = ...
            [
            -0.1;
            1
            ];
        roots_u{2,1} = ...
            [
            -0.1;
            1
            ];
        
        fxy_mtrx_exct = BuildPoly_NonSeparable(roots_f);
        gxy_mtrx_exct = BuildPoly_NonSeparable(roots_g);
        uxy_mtrx_exct = BuildPoly_NonSeparable(roots_u);
        vxy_mtrx_exct = BuildPoly_NonSeparable(roots_v);
        dxy_mtrx_exct = BuildPoly_NonSeparable(roots_d);
        
        m = 6;
        n = 7;
        t = 4;
        
        
        
    case '16'
        % same as example 13 but divided by 10
        
        roots_f{1,1} = ...
            [
            -0.2105  0   1;
            0   0   0;
            1   0   0;
            ];
        roots_f{2,1} = ...
            [
            0.2105  0   1;
            0   0   0;
            1   0   0;
            ];
        roots_f{3,1} = ...
            [
            -0.1;
            1
            ];
        roots_f{4,1} = ...
            [
            -0.1;
            1
            ];
        roots_f{5,1} = ...
            [
            0.5087652   1
            ];
        roots_f{6,1} = ...
            [
            0.5087652   1
            ];
        
        
        roots_g{1,1} = ...
            [
            -0.2105  0   1;
            0   0   0;
            1   0   0;
            ];
        roots_g{2,1} = ...
            [
            0.2105  0   1;
            0   0   0;
            1   0   0;
            ];
        roots_g{3,1} = ...
            [
            -0.2  1;
            ];
        roots_g{4,1} = ...
            [
            -0.5;
            1
            ];
        roots_g{5,1} = ...
            [
            -0.51234;
            1
            ];
        roots_g{6,1} = ...
            [
            0.5087652   1
            ];
        roots_g{7,1} = ...
            [
            0.5087652   1
            ];
        
        
        
        roots_d{1,1} = ...
            [
            -0.2105  0   1;
            0   0   0;
            1   0   0;
            ];
        roots_d{2,1} = ...
            [
            0.2105  0   1;
            0   0   0;
            1   0   0;
            ];
        roots_d{3,1} = ...
            [
            0.5087652   1
            ];
        roots_d{4,1} = ...
            [
            0.5087652   1
            ];
        
        
        
        roots_v{1,1} = ...
            [
            -0.2  1
            ];
        roots_v{2,1} = ...
            [
            -0.5;
            1
            ];
        roots_v{3,1} = ...
            [
            -0.51234;
            1
            ];
        
        roots_u{1,1} = ...
            [
            -0.1;
            1
            ];
        roots_u{2,1} = ...
            [
            -0.1;
            1
            ];
        
        fxy_mtrx_exct = BuildPoly_NonSeparable(roots_f);
        gxy_mtrx_exct = BuildPoly_NonSeparable(roots_g);
        uxy_mtrx_exct = BuildPoly_NonSeparable(roots_u);
        vxy_mtrx_exct = BuildPoly_NonSeparable(roots_v);
        dxy_mtrx_exct = BuildPoly_NonSeparable(roots_d);
        
        m = 8;
        n = 9;
        t = 6;
        
    case '17'
        % same as example 13 but divided by 10
        
        roots_f{1,1} = ...
            [
            -0.2105  0   1;
            0   0   0;
            1   0   0;
            ];
        roots_f{2,1} = ...
            [
            0.2105  0   1;
            0   0   0;
            1   0   0;
            ];
        roots_f{3,1} = ...
            [
            -0.1;
            1
            ];
        roots_f{4,1} = ...
            [
            -0.1;
            1
            ];
        roots_f{5,1} = ...
            [
            0.5087652   1
            ];
        roots_f{6,1} = ...
            [
            0.5087652   1
            ];
        roots_f{7,1} = ...
            [
            0.9728     1;
            1           0;
            ];
        
        roots_g{1,1} = ...
            [
            -0.2105  0   1;
            0   0   0;
            1   0   0;
            ];
        roots_g{2,1} = ...
            [
            0.2105  0   1;
            0   0   0;
            1   0   0;
            ];
        roots_g{3,1} = ...
            [
            -0.2  1;
            ];
        roots_g{4,1} = ...
            [
            -0.5;
            1
            ];
        roots_g{5,1} = ...
            [
            -0.51234;
            1
            ];
        roots_g{6,1} = ...
            [
            0.5087652   1
            ];
        roots_g{7,1} = ...
            [
            0.5087652   1
            ];
        roots_g{8,1} = ...
            [
            0.9728     1;
            1           0;
            ];
        
        
        
        roots_d{1,1} = ...
            [
            -0.2105  0   1;
            0   0   0;
            1   0   0;
            ];
        roots_d{2,1} = ...
            [
            0.2105  0   1;
            0   0   0;
            1   0   0;
            ];
        roots_d{3,1} = ...
            [
            0.5087652   1
            ];
        roots_d{4,1} = ...
            [
            0.5087652   1
            ];
        roots_d{5,1} = ...
            [
            0.9728     1;
            1           0;
            ];
        
        
        
        roots_v{1,1} = ...
            [
            -0.2  1
            ];
        roots_v{2,1} = ...
            [
            -0.5;
            1
            ];
        roots_v{3,1} = ...
            [
            -0.51234;
            1
            ];
        
        roots_u{1,1} = ...
            [
            -0.1;
            1
            ];
        roots_u{2,1} = ...
            [
            -0.1;
            1
            ];
        
        fxy_mtrx_exct = BuildPoly_NonSeparable(roots_f);
        gxy_mtrx_exct = BuildPoly_NonSeparable(roots_g);
        uxy_mtrx_exct = BuildPoly_NonSeparable(roots_u);
        vxy_mtrx_exct = BuildPoly_NonSeparable(roots_v);
        dxy_mtrx_exct = BuildPoly_NonSeparable(roots_d);
        
        m = 9;
        n = 10;
        t = 7;
        
        
    case '18'
        roots_f{1,1} = root_x(0.573429421);
        roots_f{2,1} = root_x(2.2175);
        roots_f{3,1} = root_x(2.2175);
        roots_f{4,1} = root_x(2.2175);
        roots_f{5,1} = root_x(0.8765);
        roots_f{6,1} = root_x(0.8765);
        roots_f{7,1} = root_x(0.1759262);
        roots_f{8,1} = root_x(1.1057853);
        
        roots_g{1,1} = root_x(0.573429421);
        roots_g{2,1} = root_x(0.7891);
        roots_g{3,1} = root_x(0.1234);
        roots_g{4,1} = root_x(0.8765);
        roots_g{5,1} = root_x(0.8765);
        roots_g{6,1} = root_x(0.0175267);
        
        roots_d{1,1} = root_x(0.573429421);
        roots_d{2,1} = root_x(0.8765);
        roots_d{3,1} = root_x(0.8765);
        
        roots_u{1,1} = root_x(2.2175);
        roots_u{2,1} = root_x(2.2175);
        roots_u{3,1} = root_x(2.2175);
        roots_u{4,1} = root_x(0.1759262);
        roots_u{5,1} = root_x(1.1057853);
        
        roots_v{1,1} = root_x(0.7891);
        roots_v{2,1} = root_x(0.1234);
        roots_v{3,1} = root_x(0.0175267);
        
        fxy_mtrx_exct = BuildPoly_NonSeparable(roots_f);
        gxy_mtrx_exct = BuildPoly_NonSeparable(roots_g);
        uxy_mtrx_exct = BuildPoly_NonSeparable(roots_u);
        vxy_mtrx_exct = BuildPoly_NonSeparable(roots_v);
        dxy_mtrx_exct = BuildPoly_NonSeparable(roots_d);
        
        m = 8;
        n = 6;
        t = 3;
        
    case '19'
        roots_f{1,1} = root_x(0.573429421);
        roots_f{2,1} = root_x(2.2175);
        roots_f{3,1} = root_x(2.2175);
        roots_f{4,1} = root_x(2.2175);
        roots_f{5,1} = root_x(0.8765);
        roots_f{6,1} = root_x(0.8765);
        roots_f{7,1} = root_x(0.1759262);
        roots_f{8,1} = root_x(1.1057853);
        % Roots in y
        roots_f{9,1} = root_x(0.424242)';
        roots_f{10,1} = root_x(0.424242)';
        roots_f{11,1} = root_x(0.424242)';
        
        roots_g{1,1} = root_x(0.573429421);
        roots_g{2,1} = root_x(0.7891);
        roots_g{3,1} = root_x(0.1234);
        roots_g{4,1} = root_x(0.8765);
        roots_g{5,1} = root_x(0.8765);
        roots_g{6,1} = root_x(0.0175267);
        % Roots in y
        roots_g{7,1} = root_x(0.424242)';
        roots_g{8,1} = root_x(0.424242)';
        roots_g{9,1} = root_x(0.424242)';
        
        
        roots_d{1,1} = root_x(0.573429421);
        roots_d{2,1} = root_x(0.8765);
        roots_d{3,1} = root_x(0.8765);
        % Roots in y
        roots_d{4,1} = root_x(0.424242)';
        roots_d{5,1} = root_x(0.424242)';
        roots_d{6,1} = root_x(0.424242)';
        
        roots_u{1,1} = root_x(2.2175);
        roots_u{2,1} = root_x(2.2175);
        roots_u{3,1} = root_x(2.2175);
        roots_u{4,1} = root_x(0.1759262);
        roots_u{5,1} = root_x(1.1057853);
        
        roots_v{1,1} = root_x(0.7891);
        roots_v{2,1} = root_x(0.1234);
        roots_v{3,1} = root_x(0.0175267);
        
        fxy_mtrx_exct = BuildPoly_NonSeparable(roots_f);
        gxy_mtrx_exct = BuildPoly_NonSeparable(roots_g);
        uxy_mtrx_exct = BuildPoly_NonSeparable(roots_u);
        vxy_mtrx_exct = BuildPoly_NonSeparable(roots_v);
        dxy_mtrx_exct = BuildPoly_NonSeparable(roots_d);
        
        m = 11;
        n = 9;
        t = 6;
        
    case '20'
        % From winkler - Two methods...
        roots_f{1,1} = root_x(-0.5161);
        roots_f{2,1} = root_x(-0.5161);
        roots_f{3,1} = root_x(-0.5161);
        roots_f{4,1} = root_x(-0.5161);
        roots_f{5,1} = root_x(-0.5161);
        
        roots_f{6,1} = root_x(-7.1052);
        roots_f{7,1} = root_x(-7.1052);
        roots_f{8,1} = root_x(-7.1052);
        roots_f{9,1} = root_x(-7.1052);
        roots_f{10,1} = root_x(-7.1052);
        
        roots_f{11,1} = root_x(-0.1132);
        roots_f{12,1} = root_x(-0.1132);
        roots_f{13,1} = root_x(-0.1132);
        
        roots_g{1,1} = root_x(-0.5161);
        roots_g{2,1} = root_x(-0.5161);
        roots_g{3,1} = root_x(-0.5161);
        roots_g{4,1} = root_x(-0.5161);
        roots_g{5,1} = root_x(-0.5161);
        
        roots_g{6,1} = root_x(-7.1052);
        roots_g{7,1} = root_x(-7.1052);
        roots_g{8,1} = root_x(-7.1052);
        roots_g{9,1} = root_x(-7.1052);
        roots_g{10,1} = root_x(-7.1052);
        
        roots_g{11,1} = root_x(2.0476);
        roots_g{12,1} = root_x(2.0476);
        roots_g{13,1} = root_x(2.0476);
        roots_g{14,1} = root_x(2.0476);
        roots_g{15,1} = root_x(2.0476);
        roots_g{16,1} = root_x(2.0476);
        roots_g{17,1} = root_x(2.0476);
        
        roots_g{18,1} = root_x(-8.8614);
        roots_g{19,1} = root_x(-8.8614);
        roots_g{20,1} = root_x(-8.8614);
        roots_g{21,1} = root_x(-8.8614);
        roots_g{22,1} = root_x(-8.8614);
        roots_g{23,1} = root_x(-8.8614);
        roots_g{24,1} = root_x(-8.8614);
        
        roots_d{1,1} = root_x(0.5161);
        roots_d{2,1} = root_x(0.5161);
        roots_d{3,1} = root_x(0.5161);
        roots_d{4,1} = root_x(0.5161);
        roots_d{5,1} = root_x(0.5161);
        
        roots_d{6,1} = root_x(7.1052);
        roots_d{7,1} = root_x(7.1052);
        roots_d{8,1} = root_x(7.1052);
        roots_d{9,1} = root_x(7.1052);
        roots_d{10,1} = root_x(7.1052);
        
        roots_u{1,1} = root_x(0.1132);
        roots_u{2,1} = root_x(0.1132);
        roots_u{3,1} = root_x(0.1132);
        
        roots_v{1,1} = root_x(2.0476);
        roots_v{2,1} = root_x(2.0476);
        roots_v{3,1} = root_x(2.0476);
        roots_v{4,1} = root_x(2.0476);
        roots_v{5,1} = root_x(2.0476);
        roots_v{6,1} = root_x(2.0476);
        roots_v{7,1} = root_x(2.0476);
        
        roots_v{8,1} = root_x(8.8614);
        roots_v{9,1} = root_x(8.8614);
        roots_v{10,1} = root_x(8.8614);
        roots_v{11,1} = root_x(8.8614);
        roots_v{12,1} = root_x(8.8614);
        roots_v{13,1} = root_x(8.8614);
        roots_v{14,1} = root_x(8.8614);
        
        m = 13;
        n = 24;
        t = 10;
        
        
        fxy_mtrx_exct = BuildPoly_NonSeparable(roots_f);
        gxy_mtrx_exct = BuildPoly_NonSeparable(roots_g);
        uxy_mtrx_exct = BuildPoly_NonSeparable(roots_u);
        vxy_mtrx_exct = BuildPoly_NonSeparable(roots_v);
        dxy_mtrx_exct = BuildPoly_NonSeparable(roots_d);
        
        
    case '21'
        % From winkler - Two methods...
        % With added roots in y
        roots_f{1,1} = root_x(-0.5161);
        roots_f{2,1} = root_x(-0.5161);
        roots_f{3,1} = root_x(-0.5161);
        roots_f{4,1} = root_x(-0.5161);
        roots_f{5,1} = root_x(-0.5161);
        
        roots_f{6,1} = root_x(-7.1052);
        roots_f{7,1} = root_x(-7.1052);
        roots_f{8,1} = root_x(-7.1052);
        roots_f{9,1} = root_x(-7.1052);
        roots_f{10,1} = root_x(-7.1052);
        
        roots_f{11,1} = root_x(-0.1132);
        roots_f{12,1} = root_x(-0.1132);
        roots_f{13,1} = root_x(-0.1132);
        
        % Roots of f in y
        roots_f{14,1} = root_x(0.7725)';
        roots_f{15,1} = root_x(0.7725)';
        roots_f{16,1} = root_x(0.7725)';
        
        % Roots of g in x
        roots_g{1,1} = root_x(-0.5161);
        roots_g{2,1} = root_x(-0.5161);
        roots_g{3,1} = root_x(-0.5161);
        roots_g{4,1} = root_x(-0.5161);
        roots_g{5,1} = root_x(-0.5161);
        
        roots_g{6,1} = root_x(-7.1052);
        roots_g{7,1} = root_x(-7.1052);
        roots_g{8,1} = root_x(-7.1052);
        roots_g{9,1} = root_x(-7.1052);
        roots_g{10,1} = root_x(-7.1052);
        
        roots_g{11,1} = root_x(2.0476);
        roots_g{12,1} = root_x(2.0476);
        roots_g{13,1} = root_x(2.0476);
        roots_g{14,1} = root_x(2.0476);
        roots_g{15,1} = root_x(2.0476);
        roots_g{16,1} = root_x(2.0476);
        roots_g{17,1} = root_x(2.0476);
        
        roots_g{18,1} = root_x(-8.8614);
        roots_g{19,1} = root_x(-8.8614);
        roots_g{20,1} = root_x(-8.8614);
        roots_g{21,1} = root_x(-8.8614);
        roots_g{22,1} = root_x(-8.8614);
        roots_g{23,1} = root_x(-8.8614);
        roots_g{24,1} = root_x(-8.8614);
        
        % Roots of g in y
        roots_g{25,1} = root_x(0.7725)';
        roots_g{26,1} = root_x(0.7725)';
        roots_g{27,1} = root_x(0.7725)';
        
        % Roots of d in x
        roots_d{1,1} = root_x(0.5161);
        roots_d{2,1} = root_x(0.5161);
        roots_d{3,1} = root_x(0.5161);
        roots_d{4,1} = root_x(0.5161);
        roots_d{5,1} = root_x(0.5161);
        
        roots_d{6,1} = root_x(7.1052);
        roots_d{7,1} = root_x(7.1052);
        roots_d{8,1} = root_x(7.1052);
        roots_d{9,1} = root_x(7.1052);
        roots_d{10,1} = root_x(7.1052);
        
        % Roots of d in y
        roots_d{11,1} = root_x(0.7725)';
        roots_d{12,1} = root_x(0.7725)';
        roots_d{13,1} = root_x(0.7725)';
        
        % Roots of u in x
        roots_u{1,1} = root_x(0.1132);
        roots_u{2,1} = root_x(0.1132);
        roots_u{3,1} = root_x(0.1132);
        
        % Roots of v in x
        roots_v{1,1} = root_x(2.0476);
        roots_v{2,1} = root_x(2.0476);
        roots_v{3,1} = root_x(2.0476);
        roots_v{4,1} = root_x(2.0476);
        roots_v{5,1} = root_x(2.0476);
        roots_v{6,1} = root_x(2.0476);
        roots_v{7,1} = root_x(2.0476);
        
        roots_v{8,1} = root_x(8.8614);
        roots_v{9,1} = root_x(8.8614);
        roots_v{10,1} = root_x(8.8614);
        roots_v{11,1} = root_x(8.8614);
        roots_v{12,1} = root_x(8.8614);
        roots_v{13,1} = root_x(8.8614);
        roots_v{14,1} = root_x(8.8614);
        
        m = 16;
        n = 27;
        t = 13;
        
        
        fxy_mtrx_exct = BuildPoly_NonSeparable(roots_f);
        gxy_mtrx_exct = BuildPoly_NonSeparable(roots_g);
        uxy_mtrx_exct = BuildPoly_NonSeparable(roots_u);
        vxy_mtrx_exct = BuildPoly_NonSeparable(roots_v);
        dxy_mtrx_exct = BuildPoly_NonSeparable(roots_d);
        
    case '22'
        % Example 2 from winkler - Two methods...
        
        % Roots of f with respect to x
        RM_mat_f_x =[...
            -9.6597 8;
            9.1433  8;
            -4.6264 9;
            10.3098 6;
            ];
        
        roots_f_x = mult_roots_x(RM_mat_f_x);
        roots_f = [roots_f_x];
        
        RM_mat_g_x = [...
            -9.6597 6;
            9.1433  6;
            0.13222 6;
            0.008550 1;
            ];
        
        roots_g_x = mult_roots_x(RM_mat_g_x);
        roots_g = [roots_g_x];
        
        RM_mat_d_x = [...
            -9.6597 6;
            9.1433  6;
        ];
    
        roots_d_x = mult_roots_x(RM_mat_d_x);
        roots_d = [roots_d_x];
        
        RM_mat_u_x = [...
            9.1433  2;
            -4.6264 9;
            10.3098 6;
            -9.6597 2;
        ];

        roots_u_x = mult_roots_x(RM_mat_u_x);
        roots_u = [roots_u_x];
        
        RM_mat_v_x = [...
            0.13222     6;
            0.008550    1;
        ];
    
        roots_v_x = mult_roots_x(RM_mat_v_x);
        roots_v = [roots_v_x];
        
        
        m = 31;
        n = 19;
        t = 12;
        
        fxy_mtrx_exct = BuildPoly_NonSeparable(roots_f);
        gxy_mtrx_exct = BuildPoly_NonSeparable(roots_g);
        uxy_mtrx_exct = BuildPoly_NonSeparable(roots_u);
        vxy_mtrx_exct = BuildPoly_NonSeparable(roots_v);
        dxy_mtrx_exct = BuildPoly_NonSeparable(roots_d);
        
    case '23'
        
        % Roots of f in terms of x
        
        % The matrix of root and multiplicity
        roots_f_x = [...
            1.2456  4;
            0.7957  3;
            0.5276  2;
            
            ];
        Arr_roots_f_x = mult_roots_x(roots_f_x);
        
        % Roots of f in terms of y
        roots_f_y = [...
            5.2476  4;
            0.1767  1;
            
            ];
        Arr_roots_f_y = mult_roots_y(roots_f_y);
        
        % Roots of f in terms of x and y
        Arr_roots_f_xy{1,1} = [...
            -2.5    0   1;
            0       0   0;
            1       0   0;
            ];
        Arr_roots_f_xy{2,1} =[
            0.75    0   1;
            0   0   0;
            1   0   0;
            ];
        roots_f = [Arr_roots_f_x; Arr_roots_f_y; Arr_roots_f_xy];
        
        fxy_mtrx_exct = BuildPoly_NonSeparable(roots_f);
        %%
        % Roots of g in terms of x
        
        % The matrix of root and multiplicity
        roots_g_x = [...
            1.2456  4;
            0.7957  3;
            0.2121  1;
            0.1212  2;
            ];
        Arr_roots_g_x = mult_roots_x(roots_g_x);
        
        % Roots of f in terms of y
        roots_g_y = [...
            5.2476  4;
            1.2765  1;
            ];
        Arr_roots_g_y = mult_roots_y(roots_g_y);
        
        % Roots of f in terms of x and y
        Arr_roots_g_xy{1,1} = [...
            -2.5    0   1;
            0       0   0;
            1       0   0;
            ];
        Arr_roots_g_xy{2,1} =[
            0.75    0   1;
            0   0   0;
            1   0   0;
            ];
        
        roots_g = [Arr_roots_g_x; Arr_roots_g_y; Arr_roots_g_xy];
        
        gxy_mtrx_exct = BuildPoly_NonSeparable(roots_g);
        %% Roots of d
        
        % The matrix of root and multiplicity
        roots_d_x = [...
            1.2456  4;
            0.7957  3;
            ];
        Arr_roots_d_x = mult_roots_x(roots_d_x);
        
        % Roots of f in terms of y
        roots_d_y = [...
            5.2476  4;
            ];
        Arr_roots_d_y = mult_roots_y(roots_d_y);
        
        % Roots of f in terms of x and y
        Arr_roots_d_xy{1,1} = [...
            -2.5    0   1;
            0       0   0;
            1       0   0;
            ]
        Arr_roots_d_xy{2,1} =[
            0.75    0   1;
            0   0   0;
            1   0   0;
            ];
        
        roots_d = [Arr_roots_d_x; Arr_roots_d_y; Arr_roots_d_xy];
        
        dxy_mtrx_exct = BuildPoly_NonSeparable(roots_d);
        
        %% Roots of u(x,y)
        
        % The matrix of root and multiplicity
        roots_u_x = [...
            0.5276  2;
            ];
        Arr_roots_u_x = mult_roots_x(roots_u_x);
        
        % Roots of f in terms of y
        roots_u_y = [...
            0.1767  1;
            ];
        Arr_roots_u_y = mult_roots_y(roots_u_y);
        
        % Roots of f in terms of x and y
        Arr_roots_u_xy{1,1} = [...
            ];
        
        roots_u = [Arr_roots_u_x; Arr_roots_u_y; Arr_roots_u_xy];
        
        uxy_mtrx_exct = BuildPoly_NonSeparable(roots_u);
        
        %% Roots of v(x,y)
        
        % The matrix of root and multiplicity
        roots_v_x = [...
            0.2121  1;
            0.1212  2;
            ];
        Arr_roots_v_x = mult_roots_x(roots_v_x);
        
        % Roots of f in terms of y
        roots_v_y = [...
            1.2765  1;
            ];
        Arr_roots_v_y = mult_roots_y(roots_v_y);
        
        % Roots of f in terms of x and y
        Arr_roots_v_xy{1,1} = [...
            ];
        
        roots_v = [Arr_roots_v_x; Arr_roots_v_y; Arr_roots_v_xy];
        
        vxy_mtrx_exct = BuildPoly_NonSeparable(roots_v);
        
        %%
        m = 21;
        n = 23;
        t = 19;
end


[r,c] = size(fxy_mtrx_exct);
m1 = r -1;
m2 = c -1;

[r,c] = size(gxy_mtrx_exct);
n1 = r -1;
n2 = c -1;

[r,c] = size(dxy_mtrx_exct);
t1 = r-1;
t2 = c-1;


end

function poly = root_x(r)
poly = [...
    -r;
    1;
    ];
end

function poly = root_y(r)
poly = [-r 1];
end


function cellArr = mult_roots_x(root_mult_mat)
% given the root and multiplicity matrix

%                    _______________
% root_mult_mat =   | root  |  mult |
%                   |       |       |
%                   |_______|_______|


count = 1;

cellArr = {}

% for each root in the array
[num_roots,~ ] = size(root_mult_mat);

for i = 1:1:num_roots
    % get the multiplicity of the root
    mult = root_mult_mat(i,2)
    % get the root
    root = root_mult_mat(i,1)
    for j = 1:1:mult
        cellArr{count,1} = root_x(root);
        count = count + 1;
    end
end

end

function cellArr = mult_roots_y(root_mult_mat)
% given the root and multiplicity matrix

%                    _______________
% root_mult_mat =   | root  |  mult |
%                   |       |       |
%                   |_______|_______|


count = 1;
cellArr = {}
% for each root in the array
[num_roots,~ ] = size(root_mult_mat);

for i = 1:1:num_roots
    % get the multiplicity of the root
    mult = root_mult_mat(i,2)
    % get the root
    root = root_mult_mat(i,1)
    for j = 1:1:mult
        cellArr{count,1} = root_y(root);
        count = count + 1;
    end
end

end