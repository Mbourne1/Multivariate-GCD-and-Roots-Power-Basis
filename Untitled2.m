
for m = 0:1:10
    for m1 = 0:1:m
        for m2 = 0:1:m
            if m1 + m2 >= m
            temp = GetNumNonZeros(m1,m2,m);
            str = sprintf('%i %i %i %i',m,m1,m2,temp);
            fprintf([str '\n']);
            end
        end
    end
end