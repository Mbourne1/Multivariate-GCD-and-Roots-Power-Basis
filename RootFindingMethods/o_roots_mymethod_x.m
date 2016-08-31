function [arr_wxy,vDegt_wx] = o_roots_mymethod_x(fxy_matrix,M)
% Given a bivariate polynomial compute the roots by differentiating with
% respect to x.

global SETTINGS

% Set the iteration number
ite = 1;

% Set the first entry of q to be the input polynomial f(x,y)
arr_fxy{1} = fxy_matrix;

% Get the total degree of f(x,y)
vDegt_arr_fxy(ite,1) = M(ite);

% Get the degree of f(x,y) with respect to x and y
[m1,m2] = GetDegree(arr_fxy{ite});
vDeg1_arr_fxy(ite,1) = m1;
vDeg2_arr_fxy(ite,1) = m2;

% Whilst the most recently calculated GCD has a degree greater than
% zero. ie is not a constant, perform a gcd calculation on it and its
% derivative.
while vDeg1_arr_fxy(ite,1) > 0
    
    if (vDeg1_arr_fxy(ite,1) == 1)
        % The derivative is a constant
        
        % The GCD is a constant
        arr_fxy{ite+1,1} = Differentiate_wrt_x(arr_fxy{ite},vDegt_arr_fxy(ite));
        
        % Deconvolve
        arr_uxy{ite+1,1} = Deconvolve_Bivariate_Single(arr_fxy{ite},arr_fxy{ite+1},vDegt_arr_fxy(ite), vDegt_arr_fxy(ite)-1);
        
        % Get total degree of d(x,y) and degree with respect to x and y
        vDegt_arr_fxy(ite+1,1) = vDegt_arr_fxy(ite)-1;
        vDeg1_arr_fxy(ite+1,1) = vDeg1_arr_fxy(ite)-1;
        vDeg2_arr_fxy(ite+1,1) = vDeg2_arr_fxy(ite);
        
        break;
    end
    
    LineBreakLarge()
    fprintf([mfilename ' : ' sprintf('GCD Calculation Loop iteration = %i \n', ite)]);
    fprintf([mfilename ' : ' sprintf('Compute GCD of f_{%i} and derivative f_{%i}\n\n',ite-1,ite-1)]);
    
    % Get the derivative of f(x,y) with respect to x.
    gxy = Differentiate_wrt_x(arr_fxy{ite},vDegt_arr_fxy(ite));
    
    % Get the total degree of f(x,y)
    m =  vDegt_arr_fxy(ite);
    
    % Get the total degree of g(x,y)
    n =  m - 1;
    
    % Get the upper and lower limit of the degree of the GCD(f, f')
    if ite > 1
        lower_lim = vDegt_arr_fxy(ite)-vNumDistinctRoots(ite-1);
        upper_lim = m-1;
    else
        lower_lim = 1;
        upper_lim = m-1;
    end
    
    fprintf([mfilename ' : ' sprintf('Minimum degree of f_{%i}: %i \n', ite, lower_lim)]);
    fprintf([mfilename ' : ' sprintf('Maximum degree of f_{%i}: %i \n', ite, upper_lim)]);
    LineBreakLarge();
    
    % GCD is only a scalar with respect to x so set equal to g(x,y).
    [arr_fxy{ite,1},~,arr_fxy{ite+1,1},arr_uxy{ite,1},arr_vxy{ite,1},t,t1,t2] = o_gcd_mymethod(arr_fxy{ite},gxy,m,n,[lower_lim, upper_lim]);
    
    % Get total degree of d(x,y) and degree with respect to x and y
    vDeg1_arr_fxy(ite+1,1) = t1;
    vDeg2_arr_fxy(ite+1,1) = t2;
    vDegt_arr_fxy(ite+1,1) = t;
    
    % Get number of distinct roots of f(ite)
    vNumDistinctRoots(ite,1) = vDegt_arr_fxy(ite) - vDegt_arr_fxy(ite+1);
    
    fprintf([mfilename ' : ' sprintf('The computed deg(GCD(f_{%i},f_{%i}) is : %i \n',ite,ite,vDegt_arr_fxy(ite+1))]);
    fprintf([mfilename ' : ' sprintf('Number of distinct roots in f_{%i} : %i \n',ite,vNumDistinctRoots(ite))]);
    fprintf([mfilename ' : ' sprintf('Degree of f_{%i} : %i \n',ite + 1, vDegt_arr_fxy(ite+1))])
    
    LineBreakLarge()
    
    % Increment the iteration number
    ite = ite + 1;
    
end

% % Obtain the series h_{i}
% Each h_{x}(i) is obtained by the deconvolution of q_{x}(i) and q_{x}(i+1)

% Get number of polynomials in the series of polynomials q_{i}
[nPolys_fxy] = size(arr_fxy,1);

% %
% %     Get h_{i}(x)
% %

% METHOD refers to method used to compute h{i}(x,y), either by
% deconvolution or from GCD triples computed above.
% From Deconvolutions
% From ux
SETTINGS.HXY_METHOD = 'From Deconvolutions';
switch SETTINGS.HXY_METHOD
    
    case 'From Deconvolutions'
        
        switch SETTINGS.DECONVOLUTION_METHOD
            
            case 'Separate' % Separate deconvolution
                
                for i = 1 : 1 : nPolys_fxy - 1 % For each pair of q_{x}(i) and q_{x}(i+1)
                    
                    % Deconvolve
                    arr_hxy{i} = Deconvolve_Bivariate_Single(arr_fxy{i}, arr_fxy{i+1},vDegt_arr_fxy(i), vDegt_arr_fxy(i+1));
                    
                end
            case 'Batch' % Batch deconvolution
                
                % Get the set of polynomials hx{i} from the deconvolution of the
                % set of polynomials fx{i}/fx{i+1}
                
                
                arr_hxy = Deconvolve_Bivariate_Batch(arr_fxy,vDegt_arr_fxy);
                
            otherwise
                error([mfilename ' : ' sprintf(' Deconvolution Method is either Separate or Batch')])
                
        end
        
    case 'From ux'
        
        arr_hxy = arr_uxy;
        
    otherwise
        error('err')
        
end

vDeg1_hx = vDeg1_arr_fxy(1:end-1) - vDeg1_arr_fxy(2:end);
vDeg2_hx = vDeg2_arr_fxy(1:end-1) - vDeg2_arr_fxy(2:end);
vDegt_hx = vDegt_arr_fxy(1:end-1) - vDegt_arr_fxy(2:end);


% %
% %
% Obtain the series of polynomials w_{x}{i}

% Each w_{x}(i) is obtained by the deconvolution of h_{x}(i) and h_{x}(i+1)

% Get number of polynomials in the array of h_{x}
[num_entries_hx] = size(arr_hxy,1);

if num_entries_hx > 1
    
    switch SETTINGS.DECONVOLUTION_METHOD
        
        case 'Separate' % Separate deconvolution
            
            for i = 1 : 1 : num_entries_hx - 1 % For each pair of q_{x}(i) and q_{x}(i+1)
                
                % Deconvolve
                arr_wxy{i,1} = Deconvolve_Bivariate_Single(arr_hxy{i},arr_hxy{i+1},vDegt_hx(i),vDegt_hx(i+1));
                
            end

        case 'Batch' % Batch deconvolution
            arr_wxy = Deconvolve_Bivariate_Batch(arr_hxy,vDegt_hx);
        otherwise
            error('err')
            
    end
    
    vDeg1_wx = vDeg1_hx(1:end-1) - vDeg1_hx(2:end);
    vDeg2_wx = vDeg2_hx(1:end-1) - vDeg2_hx(2:end);
    vDegt_wx = vDegt_hx(1:end-1) - vDegt_hx(2:end);
    
    
    % Set the final w_{x}(i+1) to be equal to h_{x}(i+1)
    arr_wxy{end+1,1} = arr_hxy{end};
    
    % Set the final degree structure
    vDeg1_wx(end) = vDeg1_hx(end);
    vDeg2_wx(end) = vDeg2_hx(end);
    vDegt_wx(end) = vDegt_hx(end);
    
else
    
    % Number of h_{i} is equal to one, so no deconvolutions to be performed
    arr_wxy{1} = arr_hxy{1};
    % Get the degree structure of h_{x,i}
    vDeg1_wx(1) = vDeg1_hx(1);
    vDeg2_wx(1) = vDeg2_hx(1);
    vDegt_wx(1) = vDegt_hx(1);
end

for i = 1:1:size(arr_wxy,1)
    fprintf([mfilename ' : ' sprintf('Roots of degree %i',i) ' \n']);
    factor = arr_wxy{i};
    
    if (length(factor) > 1)
        
        display(factor./factor(2));
        
    end
    
end

end
