function [wy,vDegt_fy] = o_roots_mymethod_y(fxy_matrix,M)
% Given a bivariate polynomial compute the roots by differentiating with
% respect to y.


% Set the iteration number
ite = 1;

% Set the first entry of q to be the input polynomial f(x,y)
fy{1} = fxy_matrix;

% Get the total degree of f(x,y)
vDegt_fy(ite) = M(ite);

% Get the degree of f(x,y) with respect to x and y
[m1,m2] = GetDegree(fy{ite});

vDeg1_fy(ite) = m1;

vDeg2_fy(ite) = m2;



% Whilst the most recently calculated GCD has a degree greater than
% zero. ie is not a constant, perform a gcd calculation on it and its
% derivative.
while vDeg2_fy(ite) > 0
    
    if (vDeg1_fy(ite) == 1)
        % The derivative is a constant
        
        % The GCD is a constant
        fy{ite+1} = Differentiate_wrt_y(fy{ite});
        uy{ite+1} = Deconvolve_Bivariate(fy{ite},fy{ite+1});
        % Get degree of d(x,y) with respect to x
        vDeg1_fy(ite+1) = 0;
        
        % Get degree of d(x,y) with respect to y
        vDeg2_fy(ite+1) = vDeg2_fy(ite);
        
        % Get total degree of d(x,y)
        vDegt_fy(ite+1) = 0;
        break;
    end
    
    LineBreakLarge()
    fprintf('GCD calculation wrt y iteration : %i \n\n',ite);
    fprintf('Compute GCD of f_{%i} and derivative f_{%i}\n\n',ite,ite);
    
    % Get the derivative of f(x,y) with respect to y.
    gxy = Differentiate_wrt_y(fy{ite});
    
    % Get the total degree of f(x,y)
    m = vDegt_fy(ite);
    
    % Get the total degree of g(x,y)
    n = m - 1;
    
    if ite > 1
        lower_lim = vDegt_fy(ite)-dy(ite-1);
        upper_lim = m-1;
    else
        lower_lim = 1;
        upper_lim = m-1;
    end
    
    fprintf([mfilename ' : ' sprintf('Minimum degree of f_{%i}: %i \n', ite+1, lower_lim)]);
    fprintf([mfilename ' : ' sprintf('Maximum degree of f_{%i}: %i \n', ite+1, upper_lim)]);
    
    % Get the GCD of f(x,y) and g(x,y)
    [fy{ite},~,fy{ite+1},uy{ite},~,t,t1,t2] = o_gcd_mymethod(fy{ite},gxy,m,n,[lower_lim,upper_lim]);
    
    % Get degree structure of q_{y}(i)
    vDegt_fy(ite+1) = t;
    
    %
    vDeg1_fy(ite+1) = t1;
    
    %
    vDeg2_fy(ite+1) = t2;
    
    % Get number of distinct roots of f(ite)
    dy(ite) = vDegt_fy(ite) - vDegt_fy(ite+1);
    
    fprintf([mfilename ' : ' sprintf('The computed deg(GCD(f_{%i},f_{%i}) is : %i \n',ite,ite,vDegt_fy(ite+1))])
    fprintf([mfilename ' : ' sprintf('Number of distinct roots in f_{%i} : %i \n',ite,dy(ite))])
    fprintf([mfilename ' : ' sprintf('Degree of f_{%i} : %i \n',ite + 1, vDegt_fy(ite+1))])
    
    LineBreakLarge()
    
    % increment the iteration number
    ite = ite + 1;
end

%% Obtain the series h_{i}

% get number of elements in the series of polynomials q_{i}
[~,num_entries_qy] = size(fy);


% %
% %     Get h_{i}(y)
% %
method = 'From Deconvolutions';
% For each pair of q_{x}(i) and q_{x}(i+1)
for i = 1 : 1 : num_entries_qy - 1
    
    % Perform the deconvolution
    switch method
        case 'From Deconvolutions'
            hy{i} = Deconvolve_Bivariate(fy{i},fy{i+1});
        case 'From ux'
            hy{i} = uy{i};
    end
    % Get the degree structure of h_{x}(i)
    deg1_hy(i) = vDeg1_fy(i) - vDeg1_fy(i+1);
    deg2_hy(i) = vDeg2_fy(i) - vDeg2_fy(i+1);
end

%% obtain the series w_{i} for the
[~,num_entries_hy] = size(hy);

if num_entries_hy > 1
    for i = 1:1:num_entries_hy-1
        wy{i}    =  Deconvolve_Bivariate(hy{i},hy{i+1});
        
        deg1_wy(i) = deg1_hy(i) - deg1_hy(i+1);
        deg2_wy(i) = deg2_hy(i) - deg2_hy(i+1);
    end
    wy{i+1} = hy{i+1};
    deg1_wy(i+1) = deg1_hy(i+1);
    deg2_wy(i+1) = deg2_hy(i+1);
else
    wy{1} = hy{1};
    deg1_wy(1) = deg1_hy(1);
    deg2_wy(1) = deg2_hy(1);
end

for i = 1:1:length(wy)
    fprintf([mfilename ' : ' sprintf('Roots of degree %i',i) ' \n']);
    factor = wy{i};
    if (length(factor) > 1)
        display(factor./factor(2));
    end
    
end

end