function [] = o_gcd(ex_num,emin,mean_method,bool_alpha_theta,low_rank_approx_method)
% o_gcd(ex_num,emin,mean_method,bool_alpha_theta,low_rank_approx_method)
%
% Calculate the GCD d(x,y) of two polynomials f(x,y) and g(x,y) taken from
% the example file.
%
% Inputs.
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



% Set the Global Variables
global PLOT_GRAPHS

SetGlobalVariables(mean_method,bool_alpha_theta,low_rank_approx_method)

EXAMPLE_TYPE = 'FromRoots';
switch EXAMPLE_TYPE
    case 'FromRoots'

% Get example polynomials
[fxy_exact, gxy_exact,...
    uxy_exact,vxy_exact,...
    dxy_exact,...
    m,m1,m2,...
    n,n1,n2,...
    t_exact,t1_exact,t2_exact] = Examples_GCD(ex_num);

   case 'FromCoefficients'
       %'-45*x*y - 15*x^3*y - 20*x*y^2 + 27*x*y^3 + 9*x^3*y^3 + 12*x*y^4';
       fxy_exact = ...
           [
           0   0   0   0  0
           0 -45 -20 +27 12
           0   0   0   0  0
           0 -15   0   9  0
           ];
       m = 6;
       m1 = 3;
       m2 = 4;
       % 45*x^2*y^2 + 15*x^2*y^3 - 27*x^2*y^4 - 9*x^2*y^5
       gxy_exact = ...
           [
           0 0  0  0   0  0
           0 0  0  0   0  0
           0 0 45 15 -27 -9
           ];
       n = 7;
       n1 = 2;
       n2 = 5;
       
       %'-70*x*y + 42*x*y^3'
       dxy_exact = ...
           [
           0 0 0 0;
           0 -70 0 42
           ]
       t_exact = 4;
       t1_exact = 1;
       t2_exact = 3;
       
end

DisplayDegreeStructure();


%% Noise

% Add noise to the coefficients of f and g
[fxy, ~] = Noise2(fxy_exact,emin);
[gxy, ~] = Noise2(gxy_exact,emin);


%% Get the GCD by zengs method
%[u,v,w] = o_gcd_zeng(fxy,gxy);


% Get GCD d(x,y) and quotient polynomials u(x,y) and v(x,y)
[uxy_calc, vxy_calc, dxy_calc,t,t1,t2] = o1(fxy,gxy,m,n);

[dxy_zeng,uxy_zeng,vxy_zeng] = o_gcd_zeng(fxy,gxy,1e-10);

%                          Results

try
    % Compare exact u(x,y) and calculated u(x,y)
    PrintCoefficients(uxy_calc,uxy_exact,'u(x,y)');
    GetDistanceBetweenPolynomials(uxy_calc,uxy_exact,'u(x,y');

    % Compare exact v(x,y) and calculated v(x,y)
    PrintCoefficients(vxy_calc,vxy_exact,'v(x,y)');
    GetDistanceBetweenPolynomials(vxy_calc,vxy_exact,'v(x,y');
catch
    fprintf('U and V are not known')
end


% Compare exact d(x,y) and calculated d(x,y)
PrintCoefficients(dxy_calc,dxy_exact,'d(x,y)');
error_dxy = GetDistanceBetweenPolynomials(dxy_calc,dxy_exact,'d(x,y');

% Compare exact d(x,y) and d(x,y) calculated by zengs method
try
PrintCoefficients(dxy_zeng,dxy_exact,'d(x,y)');
error_dxy_zeng = GetDistanceBetweenPolynomials(dxy_zeng,dxy_calc,'d(x,y)');

catch
    display(dxy_zeng)
end
PrintToFile(m,n,t,error_dxy)

% Given the two polynomials f(x,y) and g(x,y), Plot the explicit surfaces 
% z = f(x,y) and z = g(x,y).
switch PLOT_GRAPHS
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

global NOISE
global BOOL_PREPROC
global LOW_RANK_APPROXIMATION_METHOD

fullFileName = 'o_gcd_results.txt';


if exist('o_gcd_results.txt', 'file')
    fileID = fopen('o_gcd_results.txt','a');
    fprintf(fileID,'%5d \t %5d \t %5d \t %s \t %s \t %s \t %s\n',...
        m,n,t,error_dx,BOOL_PREPROC, NOISE, LOW_RANK_APPROXIMATION_METHOD);
    fclose(fileID);
else
  % File does not exist.
  warningMessage = sprintf('Warning: file does not exist:\n%s', fullFileName);
  uiwait(msgbox(warningMessage));
end




end