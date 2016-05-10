function [fxy_noisy,noise_matrix]=Noise(fxy_matrix,el,eu)
%
% Add noise to the coefficients of polynomial f
%
%
%
% Inputs
%
% fxy_matrix :
%
% el : signal to noise low limit
%
% eu : signal to noise upper limit


global SETTINGS
rng(SETTINGS.SEED)

% get the degree of input polynomial f
[m1,m2] = GetDegree(fxy_matrix);

switch nargin
    case 2 % Only one noise is specified, set upper = lower
        
        
        rp = (2*rand(m1+1,m2+1))-ones(m1+1,m2+1);
        s = rp*el;
        
        noise_matrix = fxy_matrix.*s;
        fxy_noisy = fxy_matrix + noise_matrix;
        
        
    case 3 % Specified upper and lower bound of noise
        
        y = (2*rand(m1+1,m2+1))-ones(m1+1,m2+1);
        s = eu *ones(m1+1,m2+1) -  y.*(eu-el);
        noise_matrix = fxy_matrix.*s;
        fxy_noisy = fxy_matrix + noise_matrix;
end