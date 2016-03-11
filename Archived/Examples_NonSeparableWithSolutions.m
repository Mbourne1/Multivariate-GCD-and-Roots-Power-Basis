function [fxy_mtrx_exct,gxy_mtrx_exct,...
    uxy_mtrx_exct, vxy_mtrx_exct,...
    dxy_mtrx_exct,...
    m,m1,m2,...
    n,n1,n2,...
    t,t1,t2] = Examples_NonSeparableWithSolutions(ex)

% NOTE - THESE EXAMPLES ASSUME LARGEST POWER FIRST.


switch ex
    case '0'
        % Easy Separable Example
        % Section 1.1
        % f_{x,1} and f'_{x,2}
        fxy_mtrx_exct = ...
            [
            1   3   0   -4;
            -5  -15 0   20;
            7   21  0   -28;
            -3  -9  0   12;
            ];
        
        gxy_mtrx_exct = ...
            [
            3   9   0   -12 ;
            -10 -30 0   40  ;
            7   21  0   -28 ;
            ];
        dxy_mtrx_exct = ...
            [
            1   3   0   -4  ;
            -1  -3  0   4   ;
            ];
        uxy_mtrx_exct = ...
            [
            1   ;
            -4  ;
            3   ;
            ];
        vxy_mtrx_exct = ...
            [
            3;
            -7;
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
            1   3   0   -4  ;
            -1  -3  0   4   ;
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
            [
            1;
            ];
        dxy_mtrx_exct = ...
            [
            1   3   0   -4;
            ]
        m = 4;
        n = 3;
        t = 3;
        
    case '2'
        % Section 1.2
        % Easy Separable Example
        % f_{y,1}(x,y) and f'_{y,1}(x,y) to obtain f_{y,2}
        fxy_mtrx_exct = ...
            [
            1   3   0   -4;
            -5  -15 0   20;
            7   21  0   -28;
            -3  -9  0   12;
            ];
        gxy_mtrx_exct = ...
            [
            3       6       0;
            -15     -30     0;
            21      42      0;
            -9      -18     0;
            ];
        uxy_mtrx_exct = ...
            [
            1 3 0;
            ];
        vxy_mtrx_exct = ...
            [
            3   0;
            ];
        dxy_mtrx_exct = ...
            [
            1     2 ;
            -5   -10 ;
            7    14 ;
            -3    -6 ;
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
             0    1   -5
             1   -7   10
            -1    6   -5
            ];
        gxy_mtrx_exct = ...
            [
            0    2   -10;
            1   -7    10;
            ];
        uxy_mtrx_exct = ...
            [
            0   1   ;
            1   -2  ;
            -1  1   ;
            ];
        
        vxy_mtrx_exct = ...
            [
            0   2   ;
            1   -2  ;
            ];
        
        dxy_mtrx_exct = ...
            [
            1   -5;
            ]
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
            0   0   0   1   -5;
            0   0   1   -7  10;
            0   1   -6  2   15;
            1   -7  6   28  -40;
            1   6   -1  -24 20;
            ];
        gxy_mtrx_exct = ...
            [
            0   0    0     4   -20;
            0   0    3   -21    30;
            0   2  -12     4    30;
            1   -7   6    28   -40;
            ];
        uxy_mtrx_exct = ...
            [
            0   0   0   1   ;
            0   0   1   -2  ;
            0   1   1   3   ;
            1   2   4   8   ;
            1   1   4   4   ;
            ];
        vxy_mtrx_exct = ...
            [
            0   0   0    4   ;
            0   0   3   -6   ;
            0   2   -2  -6   ;
            1   -2  -4   8   ;
            ];
        dxy_mtrx_exct = ...
            [
            1   -5;
            ];
        m = 5;
        n = 4;
        t = 1;
        
    case '5'
        % Section 3.2 - Partial wrt y
        % Non-Separable - Difficult example
        % GCD of f_{y,1} and f'_{y,1}
        fxy_mtrx_exct = ...
            [
            0   0   0   1   -5;
            0   0   1   -7  10;
            0   1   -6  2   15;
            1   -7  6   28  -40;
            1   6   -1  -24 20;
            ];
        
        gxy_mtrx_exct = ...
            [
               0       0       0        1;
               0       0       2       -7;
               0       3       -12      2;
               4       -21     12      28;
               -4      18      -2     -24;
            ];
        uxy_mtrx_exct = ...
            [
            0   0   0   1   -5  ;
            0   0   1   -6  5   ;
            0   1   -5  -4  20  ;
            1   -6  1   24  -20 ;
            ];
        vxy_mtrx_exct = ...
            [
            0   0   0   1   ;
            0   0   2   -6  ;
            0   3   -10  -4 ;
            4   -18  2   24 ;
            ];
        dxy_mtrx_exct = ...
            [
             1 ;   
            -1 ;
            ];
        
        m = 5 ;
        n = 4 ;
        
    case '6'
        % Section 4.1
        % Non-Separable - Even more difficult example
        % GCD of f_{x,1} and f'_{x,1}}
        
        fxy_mtrx_exct = ...
            [
            0        0       0       1     -13      55      -75  ;
            0        0       1     -18     120    -350      375  ;
            0        1     -17     112    -360     575      375  ;
            1      -18     121    -350     261     640     -975  ;
           -4       61    -323     579     537   -2920     2550  ;
            5      -72     346    -472    -939    3040    -2100  ;
           -2       28    -128     148     394   -1040      600  ;
            ];
        gxy_mtrx_exct = ...
            [
            0   0       0       6       -78     330     -450    ;
            0   0       0       -90     600     -1750   1875    ;
            0   4       -68     448     -1440   2300    -1500   ;
            0   -54     363     -1050   783     1920    -2925   ;
            -8  122     -646    1158    1074    -5840   5100    ;
            5   -72     346     -472    -939    3040    -2100   ;
            ];
        
        uxy_mtrx_exct = ...
            [
            0   0   0   1   ;
            0   0   1   -4  ;
            0   1   -3  1   ;
            1   -4  -2  14  ;
            -3  5   12  -20 ;
            2   -2  -8  8   ;
            ];
        vxy_mtrx_exct = ...
            [
            0   0     0     6 ;
            0   0     5   -19 ;
            0   4   -11     1 ;
            3   -11  -8    40 ;
            -5  7    20   -28 ;
            ];
        dxy_mtrx_exct = ...
            [
             1    -13   55  -75;
            -1     13  -55   75; 
            ];
        
        m = 9 ;
        n = 8 ;
        t = 4;
        
%     case '7'
%         % Section 4.2
%         fxy_mtrx_exct = ...
%             [
%             ];
%         gxy_mtrx_exct = ...
%             [
%             ];
%         uxy_mtrx_exct = ...
%             [
%             ];
%         vxy_mtrx_exct = ...
%             [
%             ];
%         dxy_mtrx_exct = ...
%             [
%             ];
%         m = ;
%         n = ;
        
    case '7'
        % Another Separable Example
        fxy_mtrx_exct = ...
            [
            1   -5  -24     -4      32   ;
            -5  25  120     20      -160 ;
            7   -35 -168    -28     224  ;
            -3  15  72      12      -96  ;
            ];
        gxy_mtrx_exct = ...
            [
            3    -15   -72    -12  96;
            -10   50   240     40 -320;
            7    -35  -168    -28  224;
            ];
        uxy_mtrx_exct = ...
            [
            1
            -4
            3
            ];
        vxy_mtrx_exct = ...
            [
            3
            -7
            ];
        dxy_mtrx_exct = ...
            [
            1   -5  -24 -4  32;
            -4  5   24  4   32;
            ];
        m = 7 ;
        n = 6 ;
        
    case '8'
        fxy_mtrx_exct = ...
            [
            0   0   1   ;
            0   2   -14 ;
            1   -14 49  ;
            ]
        gxy_mtrx_exct = ...
            [
            0   2;
            2   -14;
            ];
        
        dxy_mtrx_exct = ...
            [
            0   1;
            1   -7;
            ];
        uxy_mtrx_exct = ...
            [
            0   2;
            2   -14;
            ]
        vxy_mtrx_exct = ...
            [
            2
            ];
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