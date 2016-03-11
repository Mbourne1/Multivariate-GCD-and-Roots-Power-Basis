function [] = GivenCol_getIndex()

m = 3;
n = 1;

num_diags = (m+1) + (n+1) -1;

i_vec = 1:1:num_diags;


entr_lost_due_to_n = zeros(1,num_diags);
for i = 1:1:num_diags
    if i-(n+1) >0
        entr_lost_due_to_n(i) = i-(n+1);
    end
end

entr_lost_due_to_m = zeros(1,num_diags);
for i = 1:1:num_diags
    if i-(m+1) >0
        entr_lost_due_to_m(i) = i-(m+1);
    end
end

entr_lost_due_to_n
entr_lost_due_to_m

num_entries = i_vec - (entr_lost_due_to_n + entr_lost_due_to_m)
