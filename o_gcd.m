function [] = o_gcd(ex_num,emin,bool_preproc,low_rank_approx_method)
% o_gcd(ex_num,emin,bool_preproc,low_rank_approx_method)
%
% Calculate the GCD d(x,y) of two polynomials f(x,y) and g(x,y) taken from
% the example file.
%
%
%
% Inputs.
%
%
% ex_num  : Example Number (String)
%
% emin    : Minimum Noise level
%
% bool_preproc ('y'/'n')
%       'y' - Include Preprocessing
%       'n' - Exclude Preprocessing
%
% low_rank_approx_method ('y'/'n')
%       'Standard SNTLN'
%       'None'



% Set the Global Variables
global PLOT_GRAPHS

SetGlobalVariables(bool_preproc,low_rank_approx_method)



% Get example polynomials
[fxy_exact, gxy_exact,...
    uxy_exact,vxy_exact,...
    dxy_exact,...
    m,m1,m2,...
    n,n1,n2,...
    t_exact,t1_exact,t2_exact] = Examples_GCD(ex_num);


DisplayDegreeStructure();




%% Noise

% Add noise to the coefficients of f and g
[fxy, ~] = Noise2(fxy_exact,emin);
[gxy, ~] = Noise2(gxy_exact,emin);


%% Get the GCD by zengs method
[u,v,w] = o_gcd_zeng(fxy,gxy);


% Get GCD d(x,y) and quotient polynomials u(x,y) and v(x,y)
[uxy_calc, vxy_calc, dxy_calc,t,t1,t2] = o1(fxy,gxy,m,n);



%                          Results

% Compare exact u(x,y) and calculated u(x,y)
PrintCoefficients(uxy_calc,uxy_exact,'u(x,y)');
GetDistance(uxy_calc,uxy_exact,'u(x,y');

% Compare exact v(x,y) and calculated v(x,y)
PrintCoefficients(vxy_calc,vxy_exact,'v(x,y)');
GetDistance(vxy_calc,vxy_exact,'v(x,y');

% Compare exact d(x,y) and calculated d(x,y)
PrintCoefficients(dxy_calc,dxy_exact,'d(x,y)');
error_dxy = GetDistance(dxy_calc,dxy_exact,'d(x,y');


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

function [dist] = GetDistance(fxy_exact,fxy_computed,name)

% Normalise f(x,y) exact.
fxy_exact = NormaliseMatrix(fxy_exact);

% Normalise f(x,y) computed.
fxy_computed = NormaliseMatrix(fxy_computed);

dist = norm(fxy_exact - fxy_computed) ./ norm(fxy_exact);
fprintf('Distance of %s exact from %s computed a - b / a : \t %2.4e \n',name,name,dist);
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