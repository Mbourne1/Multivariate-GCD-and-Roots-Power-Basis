


f = ...
    [
    1 2 3;
    4 5 0;
    6 0 0;
    ]

[m1,m2] = GetDegree(f)

m = 2;

sum = 0;
prod = 1;
for i = 0:1:m
    for j = 0:1:m-i
        sum = sum + f(i+1,j+1);
        prod = prod * f(i+1,j+1);
    end
end

sum
prod

prod.^(1./nchoosek(m+2,2))

geomean(f(f~=0))