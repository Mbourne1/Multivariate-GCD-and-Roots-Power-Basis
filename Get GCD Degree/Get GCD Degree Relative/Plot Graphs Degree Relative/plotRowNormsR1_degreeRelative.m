function [] = plotRowNormsR1_degreeRelative(arr_R1_RowNorms, limits_t1, limits_t2)

lowerLimit_k1 = limits_t1(1);
upperLimit_k1 = limits_t1(2);
lowerLimit_k2 = limits_t2(1);
upperLimit_k2 = limits_t2(2);

nSubresultants_k1 = upperLimit_k1 - lowerLimit_k1 + 1;
nSubresultants_k2 = upperLimit_k2 - lowerLimit_k2 + 1;

figure_name = sprintf('Plot Row Norms of R1 from QR decomposition');
figure('name',figure_name)
hold on

for i1 = 1:1:nSubresultants_k1
    for i2 = 1:1:nSubresultants_k2
   
        k1 = lowerLimit_k1 + (i1-1);
        k2 = lowerLimit_k2 + (i2-1);
        
        vec = arr_R1_RowNorms{i1,i2};
        
        
        v_k1 = k1.* ones(length(vec),1);
        v_k2 = k2.* ones(length(vec),1);
        
        plot3(v_k1,v_k2,log10(vec));
        
    end
end

hold off

hold off





end