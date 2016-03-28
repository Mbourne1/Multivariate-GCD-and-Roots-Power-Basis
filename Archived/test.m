function [fxy_mtrx_exct] = test()

        roots_f{1,1} = ...
            [
            -0.8;
            1
            ];
        roots_f{2,1} = ...
            [
            -0.8;
            1
            ];
        roots_f{3,1} = ...
            [
            -0.3;
            1
            ];
        roots_f{4,1} = ...
            [
            -0.2    1
            ];
        roots_f{5,1} = ...
            [
            -0.2    1
            ];
        roots_f{6,1} = ...
            [
            0.1    1
            ];
              
        
        
        

fxy_mtrx_exct = BuildPoly_NonSeparable(roots_f)



end