function [] = plotSingularValues_degreeRelative(arr_SingularValues, limits_t1, limits_t2)
%
% % Inputs
%
% arr_SingularValues : Array of all Singular Values of each S_{k1,k2} 
%
% limits_t1 : upper and lower bounds on k1
%
% limits_t2 : upper and lower bounds on k2

lowerLimit_k1 = limits_t1(1);
upperLimit_k1 = limits_t1(2);
lowerLimit_k2 = limits_t2(1);
upperLimit_k2 = limits_t2(2);

nSubresultants_k1 = upperLimit_k1 - lowerLimit_k1 + 1;
nSubresultants_k2 = upperLimit_k2 - lowerLimit_k2 + 1;

figure_name = 'Plotting Singular Values';
figure('name',figure_name)
hold on

for i1 = 1:1:nSubresultants_k1
    for i2 = 1:1:nSubresultants_k2
        
        k1 = lowerLimit_k1 + (i1-1);
        k2 = lowerLimit_k2 + (i2-1);
        
        vector = arr_SingularValues{i1,i2};
        
        v_is = k1.*ones(length(vector));
        v_js = k2.*ones(length(vector));
        
        plot3(v_is,v_js,log10(vector),'*');
        
        
        
    end
end

grid on
xlabel('k_{1}');
ylabel('k_{2}');
zlabel('log10 \sigma_{k1,k2}');

hold off





end