function [] = o_gcd_batch
% Performs a batch of GCD computations for a variety of settings

%ex_num_arr = {'1'};
ex_num_arr = {'1','2','3'};
emin_arr = {1e-8};
emax_arr = {1e-12};
mean_method_arr = {'Geometric Mean Matlab Method','None'};
bool_alpha_theta_arr = {'y','n'};
low_rank_approx_method_arr = {'None','Standard STLN','Standard SNTLN'};
degree_method_arr = {'Relative','Total','Both'};

%parpool(2)

for i1 = 1:1:length(ex_num_arr)
   
    ex_num = ex_num_arr{i1};

    for i2 = 1:1:length(emin_arr)
        emin = emin_arr{i2};
        emax = emax_arr{1};
        
        for i3 = 1:1:length(mean_method_arr)
            
            mean_method = mean_method_arr{i3};
            
            for i4 = 1:1:length(bool_alpha_theta_arr)
                
                bool_alpha_theta = bool_alpha_theta_arr{i4};
                
                for i5 = 1:1:length(low_rank_approx_method_arr)
                
                    low_rank_approx_method = low_rank_approx_method_arr{i5};
                
                    parfor i6 = 1:1:length(degree_method_arr)
                        
                        degree_method = degree_method_arr{i6};
                    
                        o_gcd(ex_num,emin,emax,mean_method,bool_alpha_theta,low_rank_approx_method,degree_method)
                    
                    end
                end
            end
        end
    end

end

end