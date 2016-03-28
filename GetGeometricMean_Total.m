function [lambda] = GetGeometricMean_Total(fxy,m)

lambda = prod(fxy(fxy~=0).^ (1./nchoosek(m+2,2))) ;




end

function lambda = GetGeometricMean_Total_ByLogs(fxy,m)

% Remove zero values from f(x,y) 
fxy = abs(fxy(fxy~=0));

% Get f(x,y) in logs
fxy_log = log10(fxy);


%Multiply by nchoosek(m+2,2)
lambda_log = sum(fxy_log)

lambda = 10^(lambda_log) ^(nchoosek(m+2,2))

end