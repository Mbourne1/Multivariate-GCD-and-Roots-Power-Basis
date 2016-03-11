function [fxy_noisy,noise_matrix]=Noise(fxy_matrix,el,eu)

% Add noise to the coefficients of polynomial f

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%                           Inputs

% f :

% el : signal to noise low limit

% eu : signal to noise upper limit

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global seed
rng(seed)

% get the degree of input polynomial f
m1 = size(fxy_matrix,1) - 1;
m2 = size(fxy_matrix,2) - 1;

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