function T1 = BuildT1_Univariate(fx,n_k)
%
% % Inputs
%
% fx : Coefficient vector of polynomial f(x)
%
% n_k : Degree of polynomial v(x)


m = GetDegree_Univariate(fx);

T1 = zeros(m + n_k + 1, n_k+1);

% Each column of T_{n-k}(f(x))
for i = 0:1:n_k
    T1(i+1:m+i+1,i+1) = fx;
end

end