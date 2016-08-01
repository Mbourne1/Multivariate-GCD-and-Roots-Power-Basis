function [] = Test_Deconvolution()


% Set settings pertaining to this test

global SETTINGS
SETTINGS.PLOT_GRAPHS = 'y';
SETTINGS.MAX_ERROR_DECONVOLUTIONS = 1e-16;
SETTINGS.MAX_ITERATIONS_DECONVOLUTIONS = 50;


% Input f_{i} polynomials
x = sym('x');

% Set example number 
ex_num = '2';

switch ex_num
    case '1'
        
        % Create set of factors
        factor(1) = (x-2);
        factor(2) = (x-3);
        
        % Set multiplicity of each factor
        vMult = [7 , 12];
        
    case '2'
                                
        % Create Set of factors
        factor(1) = (x-2);
        factor(2) = (x-3.2789);
        factor(3) = (x-1.589);
        %         factor(4) = (x-0.7213);
        %         factor(5) = (x-1.5432);
        %         factor(6) = (x+5.72);
        
        % Set multiplicitiy of each factor
        vMult = [ 1 3 4 ];
end

% Get highest power of any factor
highest_pwr = max(vMult);

% Generate polynomials f_{0}(x) ,..., f_{m}(x) = 1. Where each f_{i+1}(x) is
% the f_{i+1} = GCD(f_{i},f'_{i}).
for i = 0:1:highest_pwr
    
    % Get multiplicity of each root in f_{i+1}
    mults = ((vMult - i) + abs(vMult-i)) ./2;
    
    % Get the symbolic polynomial f_{i+1}(x,y)
    f{i+1} = prod(factor.^(mults));
    
    % Get total degree of f_{i+1}(x,y)
    vDegt_fxy(i+1) = double(feval(symengine, 'degree', (f{i+1})));
end


% Get the degree structure of the polynomials h_{i}
deg_struct_h = diff(vDegt_fxy);

% Get the degree structure of the polynomials w_{i}
deg_struct_w = diff([deg_struct_h 0]);

% Get the multiplicities of the roots.
vMultiplicities = find(deg_struct_w~=0);

% Get the sequence of polynomials h_{i}(x) in symbolic form
for i = 1:1:length(f)-1
    h{i} = f{i} / f{i+1};
end

% %
% %
% Get coefficients vectors of f_{i}(x) and h_{i}(x)
for i = 1:1:length(f)
    try
        arr_fx{i,1} = sym2poly(f{i})';
        arr_hx{i,1} = sym2poly(h{i})';
    catch
        arr_fx{i,1} = 1;
    end
    
end


% %
% %
% %
Deconvolve_Bivariate_Batch_Total(arr_fx,vMultiplicities);

% %
% %
% %
Deconvolve_Bivariate_Batch_Respective(arr_fxy,vDegt_fxy);

% %
% %
% %



end