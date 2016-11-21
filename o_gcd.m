function [] = o_gcd(ex_num, emin, emax, mean_method, bool_alpha_theta, low_rank_approx_method, apf_method, degree_method)
% o_gcd(ex_num,emin,mean_method,bool_alpha_theta,low_rank_approx_method,apf_method,degree_method)
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
% apf_method
%       'None'
%       'Standard APF Linear'
%       'Standard APF Nonlinear'
%
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
% >> o_gcd('1',1e-10,1e-12, 'Geometric Mean Matlab Method', 'y','Standard STLN','Standard Nonlinear APF','Relative')


% Set the Global Variables
global SETTINGS

% add path
restoredefaultpath

addpath(...
    'APF',...
    'Build Matrices',...
    'Formatting',...
    'GCD Methods',...
    'Get Cofactors',...
    'Get Optimal Column',...
    'Get GCD Coefficients',...
    'Plotting',...
    'Results',...
    'Sylvester Matrix'...
    );
addpath(genpath('Examples'));
addpath(genpath('Get GCD Degree'));
addpath(genpath('Preprocessing'));
addpath(genpath('Low Rank Approximation'));


problem_type = 'GCD';

% % Ensure that minimum noise level is less than maximum noise level
if emin > emax
    temp = emin;
    emin = emax;
    emax = temp;
end

% Set global variables
SetGlobalVariables(problem_type, ex_num, emin, emax, mean_method, bool_alpha_theta, low_rank_approx_method, apf_method, degree_method)

% Print Parameters to screen
fprintf('INPUTS. \n')
fprintf('EXAMPLE NUMBER %s \n',ex_num)
fprintf('EMIN : %s \n',emin)
fprintf('EMAX : %s \n',emax)
fprintf('MEAN METHOD : %s \n', mean_method)
fprintf('PREPROCESSING : %s \n',bool_alpha_theta)
fprintf('LOW RANK METHOD : %s \n',low_rank_approx_method)
fprintf('APF METHOD : %s \n', apf_method)
fprintf('DEGREE METHOD : %s \n', degree_method)

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
[fxy, ~] = AddVariableNoiseToPoly(fxy_exact, emin, emax);
[gxy, ~] = AddVariableNoiseToPoly(gxy_exact, emin, emax);

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

% Get the GCD by my method
[fxy_calc, gxy_calc, dxy_calc, uxy_calc, vxy_calc,t,t1,t2] = ...
    o_gcd_mymethod(fxy,gxy,m,n,t_limits);




% Compare exact d(x,y) and calculated d(x,y)
% PrintCoefficients(dxy_calc,dxy_exact,'d(x,y)');

switch SETTINGS.DEGREE_METHOD
    case 'Relative'
        
        error.dxy = GetDistanceBetweenPolynomials(dxy_exact,dxy_calc,'d(x,y)');
        error.uxy = GetDistanceBetweenPolynomials(uxy_exact,uxy_calc,'u(x,y)');
        error.vxy = GetDistanceBetweenPolynomials(vxy_exact,vxy_calc,'v(x,y)');
        
    case 'Total'
        
        % Get d(x,y) in a matrix in terms of total degree.
        dxy_exact_total = zeros(t_exact+1,t_exact+1);
        dxy_exact_total(1:t1_exact+1,1:t2_exact+1) = dxy_exact;
        
        uxy_exact_total = zeros(m-t_exact+1, m-t_exact+1);
        uxy_exact_total(1:m1-t1_exact+1,1:m2-t2_exact+1) = uxy_exact;
        
        vxy_exact_total = zeros(n-t_exact+1, n-t_exact+1);
        vxy_exact_total(1:n1-t1_exact+1,1:n2-t2_exact+1) = vxy_exact;
        
        error.dxy = GetDistanceBetweenPolynomials(dxy_exact_total, dxy_calc, 'd(x,y)');
        error.uxy = GetDistanceBetweenPolynomials(uxy_exact_total, uxy_calc, 'u(x,y)');
        error.vxy = GetDistanceBetweenPolynomials(vxy_exact_total, vxy_calc, 'v(x,y)');
        
    case 'Both'

        error.dxy = GetDistanceBetweenPolynomials(dxy_exact,dxy_calc,'d(x,y)');
        error.uxy = GetDistanceBetweenPolynomials(uxy_exact,uxy_calc,'u(x,y)');
        error.vxy = GetDistanceBetweenPolynomials(vxy_exact,vxy_calc,'v(x,y)');
        
end



PrintToFile(m,n,t,t1,t2,error)

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





function []= PrintToFile(m,n,t,t1,t2,error)

global SETTINGS

fullFileName = sprintf('Results/Results_o_gcd%s.txt',datetime('today'));

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
       fprintf(fileID,'%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s \n',...
        datetime('now'),...
        SETTINGS.EX_NUM,...
        int2str(m),...
        int2str(n),...
        int2str(t),...
        int2str(t1),...
        int2str(t2),...
        num2str(error.uxy),...
        num2str(error.vxy),...
        num2str(error.dxy),...
        SETTINGS.MEAN_METHOD,...
        SETTINGS.BOOL_ALPHA_THETA,...
        SETTINGS.EMIN,...
        SETTINGS.EMAX,...
        SETTINGS.LOW_RANK_APPROXIMATION_METHOD,...
        int2str(SETTINGS.LOW_RANK_APPROX_REQ_ITE),...
        SETTINGS.APF_METHOD,...
        SETTINGS.DEGREE_METHOD); 
    end


    function WriteHeader()
        fprintf(fileID,'Date,EX_NUM,m,n,t,t1,t2,error_uxy,error_vxy,error_dxy,MEAN_METHOD,BOOL_ALPHA_THETA,EMIN,EMAX,LOW_RANK_APPROXIMATION_METHOD,ITERATIONS,APF_METHOD,DEGREE_METHOD\n');
    end








end