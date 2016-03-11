function [f_x_roots, f_y_roots, g_x_roots, g_y_roots,...
    u_x_roots, u_y_roots, v_x_roots, v_y_roots,...
    d_x_roots, d_y_roots, m, n, m_t, n_t,...
    t_exct, t1_exct, t2_exct] = Example_SeparableRoots(ex)



switch ex
    
    case '1'
        % Section 1.1
        % Easy Separable Example
        % GCD of f and partial wrt x
        % f_{x,1}(x,y) and f'_{x,1}(x,y) 
        
        f_x_roots = ...
            [
                1   2;
                3   1;
            ];
        f_y_roots = ...
            [
                -2  2;
                1   1;
            ];
        g_x_roots = ...
            [
                (7/3)   1;
                1   1;
                
            ];
        
        g_y_roots = ...
            [
                -2  2;
                1   1;
            ];
        
    case '2'
        
        f_x_roots = [...
                    1   1;
                    2   1;
                    3   1;
                    4   1];
                
        f_y_roots = [...
                    2   1;
                    5   1;
                    7   1;
                    ];
                
        g_x_roots = [...
                    1   1;
                    2   1;
                    4   1;
                    5   1;
                    ];
        g_y_roots = [
                    -2  1;
                    2   1;
                    12  2;
                    ];
       case '3'
        
        f_x_roots = [...
                    0.1   1;
                    0.2   1;
                    0.3   1;
                    0.4   1];
                
        f_y_roots = [...
                    0.2   1;
                    0.5   1;
                    0.7   1;
                    ];
                
        g_x_roots = [...
                    0.1   1;
                    0.2   1;
                    0.4   1;
                    0.5   1;
                    ];
        g_y_roots = [
                    -0.2  1;
                    0.2   1;
                    0.12  2;
                    ];         
    case '4'
        
        f_x_roots = [...
                    0.1   1;
                    0.2   1;
                    0.3   1;
                    0.4   3];
                
        f_y_roots = [...
                    0.2   1;
                    0.5   1;
                    0.7   1;
                    ];
                
        g_x_roots = [...
                    0.1   1;
                    0.2   1;
                    0.4   3;
                    0.5   1;
                    ];
        g_y_roots = [
                    -0.2  1;
                    0.2   1;
                    0.12  2;
                    ];   
        
    
end

    
    % Get the roots of the GCD and polynomials u and v
    [u_x_roots v_x_roots d_x_roots] = FindUVD(f_x_roots, g_x_roots);
    [u_y_roots v_y_roots d_y_roots] = FindUVD(f_y_roots, g_y_roots);
    
    % Get the roots and multiplicities of polynomial f(x)
    if ~isempty(f_x_roots)
        mult_vec_f_x = f_x_roots(:,2);
    else
        mult_vec_f_x = 0;
    end
    
    if ~isempty(f_y_roots)
        mult_vec_f_y = f_y_roots(:,2);
    else
        mult_vec_f_y = 0;
    end
    
    if ~isempty(g_x_roots)
        mult_vec_g_x = g_x_roots(:,2);
    else
        mult_vec_g_x = 0;
    end
    
    if ~isempty(g_y_roots)
        mult_vec_g_y = g_y_roots(:,2);
    else 
        mult_vec_g_y = 0;
    end
    
    if ~isempty(d_x_roots)
        mult_vec_d_x = d_x_roots(:,2);
    else
        mult_vec_d_x = 0;
    end
    
    
    if ~isempty(d_y_roots) 
        mult_vec_d_y = d_y_roots(:,2);
    else
        mult_vec_d_y = 0;
    end
    
    % Get degrees   
    m1 = sum(mult_vec_f_x);
    m2 = sum(mult_vec_f_y);
    n1 = sum(mult_vec_g_x);
    n2 = sum(mult_vec_g_y);
    t1_exct = sum(mult_vec_d_x);
    t2_exct = sum(mult_vec_d_y);
    
    % Get total degrees
    m = m1+m2;
    n = n1+n2;
    t_exct = t1_exct+t2_exct;
    m_t = m-t_exct;
    n_t = n-t_exct;
   
    

end