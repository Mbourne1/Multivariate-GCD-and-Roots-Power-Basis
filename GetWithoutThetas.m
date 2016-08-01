function fxy = GetWithoutThetas(fww,th1,th2)
% Given the coefficients of the polynomial f(w,w), divide the (i,j)-th 
% entry by theta_{1}^{i}\theta_{2}^{j}.
%
% % Inputs
%
% fww : Coefficients of f(w,w)
%
% th1 : theta_{1}
%
% th2 : theta_{2}

% Get the degree of polynomial f(x,y)
[m1,m2] = GetDegree(fww);

% Calculate f(x,y) from f(w,w)
th1_mat = diag(1./(th1.^(0:1:m1)));
th2_mat = diag(1./(th2.^(0:1:m2)));
fxy = th1_mat * fww * th2_mat;

end