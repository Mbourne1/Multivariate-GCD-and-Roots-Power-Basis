function [] = o_gcd(ex_num,emin,mean_method,bool_alpha_theta,low_rank_approx_method)
% o_gcd(ex_num,emin,mean_method,bool_alpha_theta,low_rank_approx_method)
%
% Calculate the GCD d(x,y) of two polynomials f(x,y) and g(x,y) taken from
% the example file.
%
% % Inputs.
%
%
% ex_num  : Example Number (String)
%
% emin    : Minimum Noise level
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
% % Example
% >> o_gcd('1',1e-10, 'Geometric Mean Matlab Method', 'y','None')


% Set the Global Variables
global SETTINGS

SetGlobalVariables(mean_method,bool_alpha_theta,low_rank_approx_method)


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
% Noise

% Add noise to the coefficients of f and g
[fxy, ~] = Noise2(fxy_exact,emin);
[gxy, ~] = Noise2(gxy_exact,emin);

% %
% % Get the GCD by zengs method
%[u,v,w] = o_gcd_zeng(fxy,gxy);


% Get GCD d(x,y) and quotient polynomials u(x,y) and v(x,y)
lower_limit = 1;
upper_limit = min(m,n);
t_limits = [lower_limit,upper_limit];

[uxy_calc, vxy_calc, dxy_calc,t,t1,t2] = o_gcd_mymethod(fxy,gxy,m,n,t_limits);

%[dxy_zeng,uxy_zeng,vxy_zeng] = o_gcd_zeng(fxy,gxy,1e-10);

% Results

try
    % Compare exact u(x,y) and calculated u(x,y)
    PrintCoefficients(uxy_calc,uxy_exact,'u(x,y)');
    GetDistanceBetweenPolynomials(uxy_calc,uxy_exact,'u(x,y');
    
    % Compare exact v(x,y) and calculated v(x,y)
    PrintCoefficients(vxy_calc,vxy_exact,'v(x,y)');
    GetDistanceBetweenPolynomials(vxy_calc,vxy_exact,'v(x,y');
catch
    fprintf('U and V are not known \n')
end


% Compare exact d(x,y) and calculated d(x,y)
PrintCoefficients(dxy_calc,dxy_exact,'d(x,y)');

error_dxy = GetDistanceBetweenPolynomials(dxy_calc,dxy_exact,'d(x,y');

% Compare exact d(x,y) and d(x,y) calculated by zengs method
try
    PrintCoefficients(dxy_zeng,dxy_exact,'d(x,y)');
    error_dxy_zeng = GetDistanceBetweenPolynomials(dxy_zeng,dxy_calc,'d(x,y)');
    display(dxy_zeng)
catch
    
end
PrintToFile(m,n,t,error_dxy)

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

function [] = PrintCoefficients(fxy_exact,fxy_computed,name)
fprintf('\n')
fprintf('----------------------------------------------------------------\n')
fprintf('\n')
fprintf('Comparison of %s exact and %s computed: \n',name,name)

% Normalise f(x,y) exact.
fxy_exact = NormaliseMatrix(fxy_exact);

% Normalise f(x,y) computed.
fxy_computed = NormaliseMatrix(fxy_computed);

display([GetAsVector(fxy_exact) GetAsVector(fxy_computed)]);


end




function []= PrintToFile(m,n,t,error_dx)

global SETTINGS

fullFileName = 'o_gcd_results.txt';


if exist('o_gcd_results.txt', 'file')
    fileID = fopen('o_gcd_results.txt','a');
    fprintf(fileID,'%5d \t %5d \t %5d \t %s \t %s \t %s \t %s \t %s\n',...
        m,n,t,error_dx,...
        SETTINGS.MEAN_METHOD,...
        SETTINGS.BOOL_ALPHA_THETA,...
        SETTINGS.NOISE,...
        SETTINGS.LOW_RANK_APPROXIMATION_METHOD);
    fclose(fileID);
else
    % File does not exist.
    warningMessage = sprintf('Warning: file does not exist:\n%s', fullFileName);
    uiwait(msgbox(warningMessage));
end




end