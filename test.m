
var0 = {'4','5','6'};
var1 = {'None','Standard STLN', 'Standard SNTLN'};
var2 = {'Relative','Total','Both'};

for i0 = 1:1:length(var0)
    for i1 = 1:1:length(var1)
        for i2 = 1:1:length(var2)
            
            ex_num  = var0{i0};
            low_rank_method = var1{i1};
            degree_method = var2{i2};
            
            LineBreakLarge()
            ex_num
            low_rank_method
            degree_method
            LineBreakLarge()
            LineBreakLarge()
            
            close all; clc; o_gcd(ex_num,1e-8,1e-12, 'Geometric Mean Matlab Method', 'y', low_rank_method, degree_method)
            
            
        end
    end
    
end