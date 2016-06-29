function [wx,vDegt_fx] = o_roots_mymethod_x(fxy_matrix,M)
% Given a bivariate polynomial compute the roots by differentiating with
% respect to x.


% Set the iteration number
ite = 1;

% Set the first entry of q to be the input polynomial f(x,y)
fx{1} = fxy_matrix;

% Get the total degree of f(x,y)
vDegt_fx(ite) = M(ite);

% Get the degree of f(x,y) with respect to x and y
[m1,m2] = GetDegree(fx{ite});

vDeg1_fx(ite) = m1;

vDeg2_fx(ite) = m2;

% Whilst the most recently calculated GCD has a degree greater than
% zero. ie is not a constant, perform a gcd calculation on it and its
% derivative.
while vDeg1_fx(ite) > 0
    
    if (vDeg1_fx(ite) == 1)
        % The derivative is a constant
        
        % The GCD is a constant
        fx{ite+1} = Differentiate_wrt_x(fx{ite});
        ux{ite+1} = Deconvolve_Bivariate(fx{ite},fx{ite+1});
        
        % Get degree of d(x,y) with respect to x
        vDeg1_fx(ite+1) = 0;
        
        % Get degree of d(x,y) with respect to y
        vDeg2_fx(ite+1) = vDeg2_fx(ite);
        
        % Get total degree of d(x,y)
        vDegt_fx(ite+1) = 0;
        break;
    end
    
    
    LineBreakLarge()
    fprintf([mfilename ' : ' sprintf('GCD Calculation Loop iteration = %i \n', ite)]);
    fprintf([mfilename ' : ' sprintf('Compute GCD of f_{%i} and derivative f_{%i}\n\n',ite,ite)]);
    
    % Get the derivative of f(x,y) with respect to x.
    gxy = Differentiate_wrt_x(fx{ite});
    
    % Get the total degree of f(x,y)
    m =  vDegt_fx(ite);
    
    % Get the total degree of g(x,y)
    n =  m - 1;
    
    % Get the upper and lower limit of the degree of the GCD(f, f')
    if ite > 1
        lower_lim = vDegt_fx(ite)-dx(ite-1);
        upper_lim = m-1;
    else
        lower_lim = 1;
        upper_lim = m-1;
    end
    
    fprintf([mfilename ' : ' sprintf('Minimum degree of f_{%i}: %i \n', ite+1, lower_lim)]);
    fprintf([mfilename ' : ' sprintf('Maximum degree of f_{%i}: %i \n', ite+1, upper_lim)]);
    
    % GCD is only a scalar with respect to x so set equal to g(x,y).
    [fx{ite},~,fx{ite+1},ux{ite},vx{ite},t,t1,t2] = o_gcd_mymethod(fx{ite},gxy,m,n,[lower_lim, upper_lim]);
    
    % Get degree of d(x,y) with respect to x
    vDeg1_fx(ite+1) = t1;
    
    % Get degree of d(x,y) with respect to y
    vDeg2_fx(ite+1) = t2;
    
    % Get total degree of d(x,y)
    vDegt_fx(ite+1) = t;
    
    % Get number of distinct roots of f(ite)
    dx(ite) = vDegt_fx(ite) - vDegt_fx(ite+1);
    
    fprintf([mfilename ' : ' sprintf('The computed deg(GCD(f_{%i},f_{%i}) is : %i \n',ite,ite,vDegt_fx(ite+1))]);
    fprintf([mfilename ' : ' sprintf('Number of distinct roots in f_{%i} : %i \n',ite,dx(ite))]);
    fprintf([mfilename ' : ' sprintf('Degree of f_{%i} : %i \n',ite + 1, vDegt_fx(ite+1))])
    
    LineBreakLarge()
    
    % increment the iteration number
    ite = ite + 1;
    
end




% % Obtain the series h_{i}
% Each h_{x}(i) is obtained by the deconvolution of q_{x}(i) and q_{x}(i+1)

% Get number of polynomials in the series of polynomials q_{i}
[~,num_entries_qx] = size(fx);

% %
% %     Get h_{i}(x)
% %
method = 'From Deconvolutions';
% For each pair of q_{x}(i) and q_{x}(i+1)
for i = 1 : 1 : num_entries_qx - 1
    
    % Perform the deconvolution
    switch method
        case 'From Deconvolutions'
            hx{i} = Deconvolve_Bivariate(fx{i},fx{i+1});
        case 'From ux'
            hx{i} = ux{i};
    end
    % Get the degree structure of h_{x}(i)
    deg1_hx(i) = vDeg1_fx(i) - vDeg1_fx(i+1);
    deg2_hx(i) = vDeg2_fx(i) - vDeg2_fx(i+1);
end





% % Obtain the series of polynomials w_{x}{i}
% Each w_{x}(i) is obtained by the deconvolution of h_{x}(i) and h_{x}(i+1)

% Get number of polynomials in the array of h_{x}
[~,num_entries_hx] = size(hx);

if num_entries_hx > 1
    
    % For each pair of h_{x}(i) and h_{x}(i+1) perform deconvolution
    for i = 1:1:num_entries_hx-1
        
        % Peform deconvolutions to obtain w_{x}
        
        wx{i}    = Deconvolve_Bivariate(hx{i},hx{i+1});
        
        % Get the degree structure of each w_{x}(i)
        deg1_wx(i) = deg1_hx(i) - deg1_hx(i+1);
        deg2_wx(i) = deg2_hx(i) - deg2_hx(i+1);
        
    end
    
    % Set the final w_{x}(i+1) to be equal to h_{x}(i+1)
    wx{i+1} = hx{i+1};
    
    % Set the final degree structure
    deg1_wx(i+1) = deg1_hx(i+1);
    deg2_wx(i+1) = deg2_hx(i+1);
    
else
    
    % Number of h_{i} is equal to one, so no deconvolutions to be performed
    wx{1} = hx{1};
    % Get the degree structure of h_{x,i}
    deg1_wx(1) = deg1_hx(1);
    deg2_wx(1) = deg2_hx(1);
    
end

for i = 1:1:length(wx)
    fprintf([mfilename ' : ' sprintf('Roots of degree %i',i) ' \n']);
    factor = wx{i};
    
    if (length(factor) > 1)
    
        display(factor./factor(2));
        
    end
    
end

end
