function [] = o_Deconvolution(ex_num)
% Performs a deconvolution test
%
% % Inputs
%
% ex_num : Example Number
%
% % Examples
%
% o_Deconvolution('1')

% Set settings pertaining to this test

global SETTINGS
SETTINGS.PLOT_GRAPHS = 'y';
SETTINGS.MAX_ERROR_DECONVOLUTIONS = 1e-16;
SETTINGS.MAX_ITERATIONS_DECONVOLUTIONS = 50;

% Add relevant paths
restoredefaultpath();
addpath(genpath('../Examples'));
addpath(genpath('Deconvolution'));
addpath(...
    'Build Matrices',...
    'Formatting',...
    'Preprocessing'...
    );


% Input f_{i} polynomials
x = sym('x');
y = sym('y');

% Get example
[factor_mult_arr_f] = Bivariate_Deconvolution_Examples(ex_num);

% Get vector of factors
factor = factor_mult_arr_f(:,1);

% Get vector of multiplicities of the factors
vMult = double(factor_mult_arr_f(:,2));

% Get highest power of any factor
highest_pwr = max(vMult);

% Generate polynomials f_{0}(x) ,..., f_{m}(x) = 1. Where each f_{i+1}(x) is
% the f_{i+1} = GCD(f_{i},f'_{i}).
arr_sym_fxy = cell(highest_pwr+1,1);
vDegt_arr_fxy = zeros(highest_pwr+1,1);

for i = 0:1:highest_pwr
    
    % Get multiplicity of each root in f_{i+1}
    mults = ((vMult - i) + abs(vMult-i)) ./2;
    
    % Get the symbolic polynomial f_{i+1}(x,y)
    arr_sym_fxy{i+1} = prod(factor.^(mults));
    
    % Get total degree of f_{i+1}(x,y)
    vDegt_arr_fxy(i+1) = double(feval(symengine, 'degree', (arr_sym_fxy{i+1})));
end


% Get the degree structure of the polynomials h_{i}
vDegt_arr_hxy = abs(diff(vDegt_arr_fxy));

% Get the degree structure of the polynomials w_{i}
deg_struct_w = abs(diff([vDegt_arr_hxy; 0]));

% Get the multiplicities of the roots.
vMultiplicities = find(deg_struct_w~=0);

% Get the sequence of polynomials h_{i}(x) in symbolic form
for i = 1:1:length(arr_sym_fxy)-1
    arr_sym_hxy{i,1} = arr_sym_fxy{i} / arr_sym_fxy{i+1};
end

% %
% %
% Get coefficients vectors of f_{i}(x) and h_{i}(x)
nPolys_arr_fxy = size(arr_sym_fxy,1);
nPolys_arr_hxy = size(arr_sym_hxy,1);

arr_fxy = cell(nPolys_arr_fxy,1);
arr_hxy = cell(nPolys_arr_hxy,1);


for i = 1:1:nPolys_arr_fxy
    arr_fxy{i,1} = double(rot90(coeffs(arr_sym_fxy{i},[x,y],'All'),2));
end


for i = 1:1:nPolys_arr_hxy
    arr_hxy{i,1} = double(rot90(coeffs(arr_sym_hxy{i},[x,y],'All'),2));
end

%--------------------------------------------------------------------------
% %
% %
% %
nPolys_arr_fxy = size(arr_fxy,1);

% Ensure that each f_{i}(x,y) is in a matrix of size m+1 x m+1
arr_fxy_total = cell(nPolys_arr_fxy,1);
for i = 1:1:nPolys_arr_fxy
    m = vDegt_arr_fxy(i);
    temp_mat = zeros(m+1,m+1);
    [nRows,nCols] = size(arr_fxy{i});
    temp_mat(1:nRows,1:nCols) = arr_fxy{i};
    arr_fxy_total{i} = temp_mat;
end

% Ensure that each h_{i}(x,y) is in a matrix of size (n+1,n+1) when
% comparing with the outputs of this deconvolution method.
arr_hxy_total = cell(nPolys_arr_hxy,1);
for i = 1:1:nPolys_arr_hxy
    n = vDegt_arr_hxy(i);
    temp_mat = zeros(n+1,n+1);
    [nRows,nCols] = size(arr_hxy{i});
    temp_mat(1:nRows,1:nCols) = arr_hxy{i};
    arr_hxy_total{i} = temp_mat;
end

%--------------------------------------------------------------------------
%
%   TESTING : SEPARATE - Relative Degree
%

% Initialise array
arr_hxy_Separate_Respective = cell(nPolys_arr_hxy,1);

for i = 1:1:nPolys_arr_fxy - 1
   fxy = arr_fxy{i};
   gxy = arr_fxy{i+1};
   arr_hxy_Separate_Respective{i} = Deconvolve_Bivariate_Single_Respective(fxy,gxy);
   
   
end

vErrors_Separate_Relative = GetError(arr_hxy_Separate_Respective,arr_hxy);
my_error.Separate_Relative = (norm(vErrors_Separate_Relative));
%--------------------------------------------------------------------------
%
%   TESTING : SEPARATE - Total Degree 
%

% Intialise array
arr_hxy_Separate_Total = cell(nPolys_arr_hxy,1);

for i = 1:1:nPolys_arr_fxy - 1;
   
   fxy = arr_fxy{i};
   m = vDegt_arr_fxy(i);
   gxy = arr_fxy{i+1};
   n = vDegt_arr_fxy(i+1);
   
   arr_hxy_Separate_Total{i} = Deconvolve_Bivariate_Single_Total(fxy,gxy,m,n);
end

vErrors_Separate_Total = GetError(arr_hxy_Separate_Total, arr_hxy_total);
my_error.Separate_Total = (norm(vErrors_Separate_Total));
%--------------------------------------------------------------------------
%
%   TESTING : SEPARATE - BOTH
%

% Intialise array
arr_hxy_Separate_Both = cell(nPolys_arr_hxy,1);

for i = 1:1:nPolys_arr_fxy - 1;
    
    fxy = arr_fxy{i};
    m = vDegt_arr_fxy(i);
    gxy = arr_fxy{i+1};
    n = vDegt_arr_fxy(i+1);
    
   arr_hxy_Separate_Both{i} = Deconvolve_Bivariate_Single_Both(fxy,gxy,m,n);
end

vErrors_Separate_Both = GetError(arr_hxy_Separate_Both, arr_hxy);
my_error.Separate_Both = (norm(vErrors_Separate_Both));
% -------------------------------------------------------------------------
%
%   TESTING : BATCH - TOTAL
%
arr_hxy_Batch_Total = Deconvolve_Bivariate_Batch_Total(arr_fxy_total,vDegt_arr_fxy);

% Get vector of error of each h_{i}(x,y) as a vector
vErrors_Batch_Total = GetError(arr_hxy_Batch_Total,arr_hxy_total);
my_error.Batch_Total = (norm(vErrors_Batch_Total))

%--------------------------------------------------------------------------
% 
%   TESTING : BATCH - RESPECTIVE
%
% 
arr_hxy_Batch_Relative = Deconvolve_Bivariate_Batch_Respective(arr_fxy,vDegt_arr_fxy);

% Get vector of error of each h_{i}(x,y) as a vector
vErrors_Batch_Relative = GetError(arr_hxy_Batch_Relative,arr_hxy);
my_error.Batch_Relative = (norm(vErrors_Batch_Relative));
%--------------------------------------------------------------------------
% 
%   TESTING : BATCH - BOTH  
% 

arr_hxy_Batch_Both = Deconvolve_Bivariate_Batch_Both(arr_fxy,vDegt_arr_fxy);

% Get vector of error of each h_{i}(x,y) as a vector
vErrors_Batch_Both = GetError(arr_hxy_Batch_Both, arr_hxy);
my_error.Batch_Both = (norm(vErrors_Batch_Both));

%-------------------------------------------------------------------------

PrintToResultsFile(ex_num,bool_preproc,emin,my_error);


end

function [] = PrintToResultsFile(ex_num, bool_preproc, emin, my_error)
%
% % Inputs
%
% ex_num :Example Number
%
% bool_preproc


fullFileName = sprintf('Deconvolution/Results/Results_o_deconvolutions%s.txt',datetime('today'));

% If file already exists append a line
if exist(fullFileName, 'file')
    
    fileID = fopen(fullFileName,'a');
    WriteNewLine()
    fclose(fileID);
    
else % File doesnt exist so create it
    
    fileID = fopen( fullFileName, 'wt' );
    WriteHeader()
    WriteNewLine()
    fclose(fileID);
    
end

    function WriteNewLine()
        
        %
        fprintf(fileID,'%s,%s,%s,%s,%s,%s,%s,%s,%s,%s \n',...
            datetime('now'),...
            ex_num,...
            bool_preproc,...
            emin,...
            my_error.Separate_Total,...
            my_error.Separate_Relative,...
            my_error.Separate_Both ,...
            my_error.Batch_Total,...
            my_error.Batch_Relative,...
            my_error.Batch_Both...
            );
        
    end

    function WriteHeader()
        
        fprintf(fileID,'DATE,EX_NUM,BOOL_PREPROC,NOISE,Sep_Total,Sep_Relative,Sep_Both,Batch_Total,Batch_Relative,Batch_Both \n');
        
    end


end


function v_errors =  GetError(arr_hxy_comp,arr_hxy)

% Get number of polynomials h_{i}(x,y) in array
nPolys_arr_hxy = size(arr_hxy_comp,1);

for i = 1:1:nPolys_arr_hxy
    
    hxy_comp_norm = arr_hxy_comp{i} ./ arr_hxy_comp{i}(1,1);
    hxy_exact = arr_hxy{i}./ arr_hxy{i}(1,1);
    
    err = abs(hxy_exact - hxy_comp_norm) ./ hxy_exact;
    
    % remove nan values
    err(isnan(err)) = [];
    
    % remove inf values
    err(isinf(err)) = [];
    
    v_errors(i) = norm(err);
   
end





end