function [arr_wxy,vDegt_wy] = o_roots_mymethod_y(fxy_matrix,M)
% Given a bivariate polynomial compute the roots by differentiating with
% respect to y.

global SETTINGS

% Set the iteration number
ite = 1;

% Set the first entry of q to be the input polynomial f(x,y)
arr_fxy{1,1} = fxy_matrix;

% Get the total degree of f(x,y)
vDegt_arr_fxy(ite) = M(ite);

% Get the degree of f(x,y) with respect to x and y
[m1,m2] = GetDegree(arr_fxy{ite});
vDeg1_arr_fxy(ite) = m1;
vDeg2_arr_fxy(ite) = m2;

% Whilst the most recently calculated GCD has a degree greater than
% zero. ie is not a constant, perform a gcd calculation on it and its
% derivative.
while vDeg2_arr_fxy(ite) > 0
    
    if (vDeg2_arr_fxy(ite) == 1)
        % The derivative is a constant
        
        % The GCD is a constant
        arr_fxy{ite+1,1} = Differentiate_wrt_y(arr_fxy{ite},vDegt_arr_fxy(ite));
                
        % Deconvolve
        arr_uxy{ite+1,1} = Deconvolve_Bivariate_Single(arr_fxy{ite},arr_fxy{ite+1},vDegt_arr_fxy(ite),vDegt_arr_fxy(ite)-1);
        
        % Get total degree of d(x,y) and degree with respect to x and y
        vDeg1_arr_fxy(ite+1) = vDeg1_arr_fxy(ite)-1;
        vDeg2_arr_fxy(ite+1) = vDeg2_arr_fxy(ite);
        vDegt_arr_fxy(ite+1) = vDegt_arr_fxy(ite)-1;
        break;
    end
    
    LineBreakLarge()
    fprintf([mfilename ' : ' sprintf('GCD calculation wrt y iteration : %i \n\n',ite)]);
    fprintf([mfilename ' : ' sprintf('Compute GCD of f_{%i} and derivative f_{%i}\n\n',ite,ite)]);
    
    % Get the derivative of f(x,y) with respect to y.
    gxy = Differentiate_wrt_y(arr_fxy{ite},vDegt_arr_fxy(ite));
    
    % Get the total degree of f(x,y)
    m = vDegt_arr_fxy(ite);
    
    % Get the total degree of g(x,y)
    n = m - 1;
    
    if ite > 1
        lower_lim = vDegt_arr_fxy(ite)-dy(ite-1);
        upper_lim = m-1;
    else
        lower_lim = 1;
        upper_lim = m-1;
    end
    
    
    fprintf([mfilename ' : ' sprintf('Minimum degree of f_{%i}: %i \n', ite+1, lower_lim)]);
    fprintf([mfilename ' : ' sprintf('Maximum degree of f_{%i}: %i \n', ite+1, upper_lim)]);
    
    % Get the GCD of f(x,y) and g(x,y)
    [arr_fxy{ite,1},~,arr_fxy{ite+1,1},arr_uxy{ite,1},~,t,t1,t2] = o_gcd_mymethod(arr_fxy{ite},gxy,m,n,[lower_lim,upper_lim]);
    
    % Get total degree of d(x,y) and degree with respect to x and y
    vDeg1_arr_fxy(ite+1,1) = t1;
    vDeg2_arr_fxy(ite+1,1) = t2;
    vDegt_arr_fxy(ite+1) = t;
    
    % Get number of distinct roots of f(ite)
    dy(ite,1) = vDegt_arr_fxy(ite) - vDegt_arr_fxy(ite+1);
    
    fprintf([mfilename ' : ' sprintf('The computed deg(GCD(f_{%i},f_{%i}) is : %i \n',ite,ite,vDegt_arr_fxy(ite+1))])
    fprintf([mfilename ' : ' sprintf('Number of distinct roots in f_{%i} : %i \n',ite,dy(ite))])
    fprintf([mfilename ' : ' sprintf('Degree of f_{%i} : %i \n',ite + 1, vDegt_arr_fxy(ite+1))])
    
    LineBreakLarge()
    
    % Increment the iteration number
    ite = ite + 1;
end

% %
% %
% % Obtain the series h_{i}


% Get number of elements in the series of polynomials q_{i}
[nPolys_arr_fxy] = size(arr_fxy,1);

% %
% %     Get h_{i}(y)
% %

% METHOD refers to method used to compute h{i}(x,y), either by
% deconvolution or from GCD triples computed above.
% From Deconvolutions
% From uy
SETTINGS.HXY_METHOD = 'From Deconvolutions';
switch SETTINGS.HXY_METHOD
    
    case 'From Deconvolutions'
        
        switch SETTINGS.DECONVOLUTION_METHOD
            
            case 'Separate'
                % For each pair of q_{x}(i) and q_{x}(i+1)
                for i = 1 : 1 : nPolys_arr_fxy - 1
                    
                    % Perform deconvolutions
                    arr_hxy{i,1} = Deconvolve_Bivariate_Single(arr_fxy{i},arr_fxy{i+1},vDegt_arr_fxy(i),vDegt_arr_fxy(i+1));
                    
                end
            case 'Batch'
                
                % Get the set of polynomials hy{i} from the deconvolution of the
                % set of polynomials fy{i}/fy{i+1}
                arr_hxy = Deconvolve_Bivariate_Batch(arr_fxy,vDegt_arr_fxy);
                 
            otherwise
                error('err')
                
        end
        
    case 'From ux'
        
        arr_hxy = arr_uxy;
        
    otherwise
        error('err');
        
end

vDeg1_hy = vDeg1_arr_fxy(1:end-1) - vDeg1_arr_fxy(2:end);
vDeg2_hy = vDeg2_arr_fxy(1:end-1) - vDeg2_arr_fxy(2:end);
vDegt_hy = vDegt_arr_fxy(1:end-1) - vDegt_arr_fxy(2:end);


% %
% %
% Obtain the series of polynomials w_{y}{i} for the

% Each w_{y}(i) is obtained by the deconvolution of h_{y}(i) and h_{y}(i+1)

% Get number of polynomials in the array of h_{y}
[nPolys_hxy] = size(arr_hxy,1);

if nPolys_hxy > 1
    
    switch SETTINGS.DECONVOLUTION_METHOD
        
        case 'Separate' % Separate deconvolution
            
            for i = 1 : 1 : nPolys_hxy - 1 % For each pair of h_{y}(i) and h_{y}(i+1)
                
                % Deconvolve
                arr_wxy{i,1} = Deconvolve_Bivariate_Single(arr_hxy{i},arr_hxy{i+1},vDegt_hy(i),vDegt_hy(i+1));
                
            end

        case 'Batch' % Batch deconvolution
            arr_wxy = Deconvolve_Bivariate_Batch(arr_hxy,vDegt_hy);
        otherwise
            error('err')
            
    end
    
    vDeg1_wy = vDeg1_hy(1:end-1) - vDeg1_hy(2:end);
    vDeg2_wy = vDeg2_hy(1:end-1) - vDeg2_hy(2:end);
    vDegt_wy = vDegt_hy(1:end-1) - vDegt_hy(2:end);
    
    
    % Set the final w_{x}(i+1) to be equal to h_{x}(i+1)
    arr_wxy{end+1,1} = arr_hxy{end};
    
    % Set the final degree structure
    vDeg1_wy(end) = vDeg1_hy(end);
    vDeg2_wy(end) = vDeg2_hy(end);
    vDegt_wy(end) = vDegt_hy(end);
    
else
    
    % Number of h_{i} is equal to one, so no deconvolutions to be performed
    arr_wxy{1} = arr_hxy{1};
    % Get the degree structure of h_{x,i}
    vDeg1_wy(1) = vDeg1_hy(1);
    vDeg2_wy(1) = vDeg2_hy(1);
    vDegt_wy(1) = vDegt_hy(1);
end

for i = 1:1:length(arr_wxy)
    fprintf([mfilename ' : ' sprintf('Roots of degree %i',i) ' \n']);
    factor = arr_wxy{i};
    if (length(factor) > 1)
        display(factor./factor(2));
    end
    
end

end