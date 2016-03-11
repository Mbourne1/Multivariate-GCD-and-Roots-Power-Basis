function [f]=BuildPoly_Pwr(A)
% Obtain polynomial coefficients in the power form, given a set of roots 
% and multiplicities. Coefficients are given in order of descending power.
%
%       Inputs
%
%   A 

% Calculate the number of distinct roots of the polynomial.
r = size(A,1);

% Convolve each factor, which is defined by a row of A, separately.
% A(k,1) stores the value of the root, and A(k,2) stores its multiplicity.

f = 1;
% for each unique root 1,...,r
for k = 1:1:r
    
    w = pwr_conv(A(k,1),A(k,2));
    f = conv(f,w) ;
end

% Obtain coefficients so that the leading coefficient has the highest
% power.
f = fliplr(f);

f = f';
end

function [poly] = pwr_conv(root,mult)
% This function convolves the vector [-r 1] with itself m times, and
% returns the corresponding polynomial.
% 
%
%                           Inputs:
%
%
% root :   Root
%
% mult :   Multiplicity of root
%
%
%                           Outputs:
%
%
% t :   Vector which stores the result from this convolution.
%



% If the multiplicity of the root is only 1, then output the polynomial 
% [ -r , 1]
if mult==1
    poly = [-root,1];
else
    % Perform polynomial multiplication m times where m is the multiplicity
    % of the root.
    
    % Build polynomial of the simple root
    q = [-root,1];
    
    % Build the starting polynomial 
    poly = [-root,1];
    
    for k = 2 : 1 : mult
        % Multiply the polynomial by the factor. (x-r)
        poly = conv(poly,q);
    end
end

end
