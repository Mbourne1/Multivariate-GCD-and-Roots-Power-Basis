function [partial_fxy] = Differentiate_wrt_y(fxy_matrix)
% Given the polynomial fxy in matrix form
%
%           1   y   y^{2}
%          _           _
%       1 |             |
%       x |             |
%   x^{2} |_           _|
%
% Obtain the partial derivative with respect to y
%

% Get the degree of fxy with respect to x and y
[m1,m2] = GetDegree(fxy_matrix);

% create a vector to multiply coefficients by

mult_vec = (0:1:m2);

% multiply the rows of fxy by the multiplication vector
partial_fxy =  fxy_matrix * diag(mult_vec); 

partial_fxy(:,1) = [];

end