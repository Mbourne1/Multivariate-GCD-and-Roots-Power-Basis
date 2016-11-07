function [] = o_gcd(ex_num,emin,emax,mean_method,bool_alpha_theta,low_rank_approx_method,degree_method)
% o_gcd(ex_num,emin,mean_method,bool_alpha_theta,low_rank_approx_method)
%
% Calculate the GCD d(x,y) of two polynomials f(x,y) and g(x,y) taken from
% the example file.
%
% % Inputs.
%
% ex_num  : Example Number (String)
%
% emin    : Minimum Noise level
%
% emax : Maximum signal to noise ratio
%
% mean_method :
%       'Geometric Mean Matlab Method'
%       'None'
%
% bool_alpha_theta ('y'/'n')
%       'y' - Include Preprocessing
%       'n' - Exclude Preprocessing
%
% low_rank_approx_method ('y'/'n')
%       'Standard SNTLN'
%       'None'
%
% degree_method
%       'Relative' : Define polynomials in terms of degree with respect to
%                    x and y, so matrices of coefficients are rectangular.
%       'Total' :   Define polynomials in terms of total degree, so matrices
%                   of coefficients are square, where lower right triangle
%                   are all zero.
%
%       'Both' :    Combination of both of the above, and typically gives best
%                   results.
%
% % Example
% >> o_gcd('1',1e-10,1e-12, 'Geometric Mean Matlab Method', 'y','Standard STLN','Relative')


% Set the Global Variables
global SETTINGS

% add path
restoredefaultpath

addpath(...
    'Build Matrices',...
    'Examples',...
    'Examples/Examples_GCD',...
    'Formatting',...
    'GCD Methods',...
    'Get Cofactors',...
    'Get Optimal Column',...
    'Get GCD Coefficients',...
    'Plotting',...
    'Sylvester Matrix'...
    );
addpath(genpath('Get GCD Degree'));
addpath(genpath('Preprocess'));
addpath(genpath('Low Rank Approximation'));


problem_type = 'GCD';

SetGlobalVariables(problem_type,ex_num,emin,emax,mean_method,bool_alpha_theta,low_rank_approx_method,degree_method)


% Get example polynomials
[fxy_exact, gxy_exact,...
    uxy_exact,vxy_exact,...
    dxy_exact,...
    m,m1,m2,...
    n,n1,n2,...
    t_exact,t1_exact,t2_exact] = Examples_GCD(ex_num);




DisplayDegreeStructure();

% %
% %
% Add Noise

% Add noise to the coefficients of f and g
[fxy, ~] = AddNoiseToPoly2(fxy_exact,emin);
[gxy, ~] = AddNoiseToPoly2(gxy_exact,emin);

% %
% % Get the GCD by zengs method
%[u,v,w] = o_gcd_zeng(fxy,gxy);


% Get GCD d(x,y) and quotient polynomials u(x,y) and v(x,y)
lower_limit = 1;
upper_limit = min(m,n);
t_limits = [lower_limit,upper_limit];



switch SETTINGS.DEGREE_METHOD
    case 'Total'
        fxy_matrix_padd = zeros(m+1,m+1);
        gxy_matrix_padd = zeros(n+1,n+1);

        [r,c] = size(fxy);
        fxy_matrix_padd(1:r,1:c) = fxy;
        
        [r,c] = size(gxy);
        gxy_matrix_padd(1:r,1:c) = gxy;
        
        fxy = fxy_matrix_padd;
        gxy = gxy_matrix_padd;
        
    case 'Relative'
        
    case 'Both'
end


[fxy_calc, gxy_calc, dxy_calc, uxy_calc, vxy_calc,t,t1,t2] = o_gcd_mymethod(fxy,gxy,m,n,t_limits);




% Compare exact d(x,y) and calculated d(x,y)
% PrintCoefficients(dxy_calc,dxy_exact,'d(x,y)');
try
    switch SETTINGS.DEGREE_METHOD
        case 'Relative'
            
            error_dxy = GetDistanceBetweenPolynomials(dxy_exact,dxy_calc,'d(x,y');
            
        case 'Total'
           
            % Get d(x,y) in a matrix in terms of total degree.
            dxy_exact_total = zeros(t_exact+1,t_exact+1);
            dxy_exact_total(1:t1_exact+1,1:t2_exact+1) = dxy_exact;
            error_dxy = GetDistanceBetweenPolynomials(dxy_exact_total,dxy_calc,'d(x,y');
            
        case 'Both'
            
            error_dxy = GetDistanceBetweenPolynomials(dxy_exact,dxy_calc,'d(x,y');
    end
    
catch err
    fprintf([mfilename ' : ' err.message '\n']);
    error_dxy = 9999999;
end

PrintToFile(m,n,t,t1,t2,error_dxy)

% Given the two polynomials f(x,y) and g(x,y), Plot the explicit surfaces
% z = f(x,y) and z = g(x,y).
switch SETTINGS.PLOT_GRAPHS
    case 'y'
        
        surfaces{1} = fxy;
        surfaces{2} = gxy;
        surfaces{3} = dxy_calc;
        
        PlotExplicitSurfaces(surfaces);
        
    case 'n'
    otherwise
        error('error - plotgraphs is either y or n');
end




end





function []= PrintToFile(m,n,t,t1,t2,error_dx)

global SETTINGS

fullFileName = 'Results_o_gcd.txt';


if exist('Results_o_gcd.txt', 'file')
    fileID = fopen('Results_o_gcd.txt','a');
    fprintf(fileID,'%s, \t %s, \t %s, \t %s, \t %s, \t %s, \t %s, \t %s, \t %s, \t %s, \t %s, \t %s, \t %s, \t %s\n',...
        datetime('now'),...
        SETTINGS.EX_NUM,...
        int2str(m),...
        int2str(n),...
        int2str(t),...
        int2str(t1),...
        int2str(t2),...
        error_dx,...
        SETTINGS.MEAN_METHOD,...
        SETTINGS.BOOL_ALPHA_THETA,...
        SETTINGS.EMIN,...
        SETTINGS.EMAX,...
        SETTINGS.LOW_RANK_APPROXIMATION_METHOD,...
        SETTINGS.DEGREE_METHOD);
    fclose(fileID);
else
    % File does not exist.
    warningMessage = sprintf('Warning: file does not exist:\n%s', fullFileName);
    uiwait(msgbox(warningMessage));
end















end