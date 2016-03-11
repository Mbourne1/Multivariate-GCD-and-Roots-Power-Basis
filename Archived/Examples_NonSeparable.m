function [fxy_mtrx_exct,gxy_mtrx_exct] = Examples_NonSeparable(ex)

switch ex
    case '0'
        % Easy Separable Example
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
        
        
% % Bivariate NonSeparable Example - Section 2

    case '3'
        % Section 2.1
        %f_{x,1} and f'_{x,1}
        fxy_mtrx_exct = ...
            [
                0   1   -5
                1   -7  10
                -1  6   -5
            ];
        gxy_mtrx_exct = ...
            [
                0 2 -10;
                1 -7 10;
            ];

% % Bivariate NonSeparable - Difficult Example - Section 3       
        
    case '4'
        % Bivariate Case
        % Section 3.1 
        % f_{x,1} and f'_{x,1}
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
            0   0   0     4   -20;
            0   0   3   -21    30;
            0   0   0   -12     4;
            1   -7  6    28   -40;
            ];
        
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
            0   0       0       0        1;
            0   0       0       2       -7;
            0   0       3       -12      2;
            0   4       -21     12      28;
            0   -4      18      -2     -24;
            ];
    case '6'
        % Section 4.1
        % Non-Separable - Even more difficult example
        % GCD of f_{x,1} and f'_{x,1}}
        
        fxy_mtrx_exct = ...
            [
            0       0       0       1       -13      55      -75    ;
            0       0       1       -18     120     -350    375     ;
            0       1       -17     112     -360    575     375     ;
            1       -18     121     -350    261     640     -975    ;
            -4      61      -323    579     537     -2920   2550    ;
            5       -72     346     -472    -939    3040    -2100   ;
            -2      28      -128    148     394     -1040   600;
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
end
end