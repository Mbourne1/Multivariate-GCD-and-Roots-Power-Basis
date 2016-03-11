function [u_x_roots, v_x_roots, d_x_roots] = FindUVD(roots_f_x, g_x_roots)
% Given the separable roots of f and g. Return the roots of u v and d.

% initialise empty vectors for d, u and v

d_x_roots = [];
u_x_roots = [];
v_x_roots = [];

% for each root in f
for i = 1:1:size(roots_f_x)
    root = roots_f_x(i,1);
    % Get the multiplicity of the root in f
    mult_f = roots_f_x(i,2);
    % find if the root exists in g
    [r , ~] = find(g_x_roots(:,1) == root);
    
    % if the root of f doesnt exist in g
    if isempty(r)
        % add the route to the polynomial u
        rm = [root mult_f];
        u_x_roots = [u_x_roots; rm];
    else
        
        
        % Given the row of the root as appearing in g, find its
        % multiplicity in g
        mult_g = g_x_roots(r,2);
        
        % If the multiplicities are the same, add to the gcd vector
        if (mult_f == mult_g)
            
            % Assign mult_d = mult_f = mult_g
            mult_d = mult_f;
            
            % Create vector of root and mult
            rm = [root mult_d];
            
            % Add to the matrix of roots for the GCD
            d_x_roots = [d_x_roots ; rm];
            
            % Remove the root from g
            g_x_roots(r,:) = [];
            
            
        elseif (mult_f > mult_g)
            
            % Assign mult_d to the lower of mult_f or mult_g
            mult_d = mult_g;
            
            % Create Vector of root and mult for insertion to GCD
            rm = [root mult_d];
            
            % Add to the matrix of roots for the GCD
            d_x_roots = [d_x_roots; rm ];
            
            % % Add the roots to u matrix
            % get the multiplicity of the root in u
            mult_u = mult_f - mult_d;
            
            % Create vector of root and mult for insertion to u
            rm = [root mult_u];
            
            % Add to the matrix of roots for the quotient u
            u_x_roots = [u_x_roots; rm];
            
            % Remove the root from g
            g_x_roots(r,:) = [];
            
            
        elseif (mult_f < mult_g)
            
            % assign mult_d
            mult_d = mult_f;
            
            % Create vector of root and mult for insertion to GCD
            rm = [root mult_d];
            
            % Add to the matrix of roots for the GCD
            d_x_roots = [d_x_roots; rm];
            
            % Get the multiplicity of the root in v
            mult_v = mult_g - mult_d;
            
            % Create vector of root and mult for insertion to v
            rm = [roots_v_x; rm];
            
            % Remove the root from g
            g_x_roots(r,:) = [];
        end
    end
    
end

    % Any remaining roots in g should go immediately to v, since they do
    % not exist in f.
    v_x_roots = [v_x_roots ; g_x_roots];
    