function [fxy_matrix_Pwr, gxy_matrix_Pwr,...
    uxy_matrix_Pwr,vxy_matrix_Pwr,dxy_matrix_Pwr...
    m,n,n_t,m_t,...
    t_exct,t1_exct,t2_exct] = GetCoeff(ex)
% Given the roots of the polynomials f(x,y) g(x,y), u(x,y), v(x,y) and
% d(x,y) return the matrices of their coefficients. The entries of these
% matrices are of the form a_{i,j} x^{m1-i} y^{m2-j}. The first entry has 
% the maximum powers x^{m1}y^{m2}, and the last bottom right entry has 
% powers x^{0}y^{0}.

% Get the polynomial roots and multiplicities
[roots_f_x, roots_f_y, ...
    roots_g_x, roots_g_y, ...
    roots_u_x, roots_u_y, ...
    roots_v_x, roots_v_y, ...
    roots_d_x, roots_d_y,...
    m, n, m_t, n_t,...
    t_exct,t1_exct,t2_exct] = Example_SeparableRoots(ex);

% % BUILD THE POLYNOMIALS
% Get coefficients of polynomial f 
f_x_poly_pwr = BuildPoly_Pwr(roots_f_x);
f_y_poly_pwr = BuildPoly_Pwr(roots_f_y);
fxy_matrix_Pwr = f_x_poly_pwr * f_y_poly_pwr';

% get the degrees m1 and m2
m1 = length(f_x_poly_pwr) - 1;
m2 = length(f_y_poly_pwr) - 1;

fxy_vec_Pwr = [];
for i = 0:1:(m1+m2)
    for j = 0:1:i
        
        row = m1 -(m1-j)+1;
        col = m2 -(m2-(i-j)) + 1;
        
        % if within the number of columns and number of rows
        if col <= (m2 + 1) && row <= (m1 + 1)
            temp = fxy_matrix_Pwr(m1-(m1-j)+1,m2-(m2 - (i-j))+1);
            fxy_vec_Pwr = [fxy_vec_Pwr ; temp];
        end
        
    end
end

% Get coefficients of polynomial g
g_x_poly_pwr = BuildPoly_Pwr(roots_g_x);
g_y_poly_pwr = BuildPoly_Pwr(roots_g_y);
gxy_matrix_Pwr = g_x_poly_pwr * g_y_poly_pwr';

n1 = length(g_x_poly_pwr)-1;
n2 = length(g_y_poly_pwr)-1;

gxy_vec_Pwr = [];
for i = 0:1:(n1+n2)
    for j = 0:1:i
        
        row = n1 -(n1-j)+1;
        col = n2 -(n2-(i-j)) + 1;
        
        % if within the number of columns and number of rows
        if col <= (n2 + 1) && row <= (n1 + 1)
            temp = gxy_matrix_Pwr(n1-(n1-j)+1,n2-(n2 - (i-j))+1);
            gxy_vec_Pwr = [gxy_vec_Pwr ; temp];
        end
        
    end
end

% Get coefficients of polynomial u
u_x_poly_pwr = BuildPoly_Pwr(roots_u_x);
u_y_poly_pwr = BuildPoly_Pwr(roots_u_y);
uxy_matrix_Pwr = u_x_poly_pwr * u_y_poly_pwr';

m1_t1 = length(u_x_poly_pwr)-1;
m2_t2 = length(u_y_poly_pwr)-1;

uxy_vec_Pwr = [];
for i = 0:1:(m1_t1+m2_t2)
    for j = 0:1:i
        
        row = m1 -(m1_t1-j)+1;
        col = m2 -(m2_t2-(i-j)) + 1;
        
        % if within the number of columns and number of rows
        if col <= (m2_t2 + 1) && row <= (m1_t1 + 1)
            temp = fxy_matrix_Pwr(m1_t1-(m1_t1-j)+1,m2_t2-(m2_t2 - (i-j))+1);
            uxy_vec_Pwr = [uxy_vec_Pwr ; temp];
        end
    end
end

% Get coefficients of polynomial u
v_x_poly_pwr = BuildPoly_Pwr(roots_v_x);
v_y_poly_pwr = BuildPoly_Pwr(roots_v_y);
vxy_matrix_Pwr = v_x_poly_pwr * v_y_poly_pwr';

n1_t1 = length(v_x_poly_pwr)-1;
n2_t2 = length(v_y_poly_pwr)-1;

vxy_vec_Pwr = [];
for i = 0:1:(n1_t1+n2_t2)
    for j = 0:1:i
  
        row = n1_t1 -(n1_t1-j)+1;
        col = n2_t2 -(n2_t2-(i-j)) + 1;
        
        % if within the number of columns and number of rows
        if col <= (n2_t2 + 1) && row <= (n1_t1 + 1)
            temp = fxy_matrix_Pwr(n1_t1-(n1_t1-j)+1,n2_t2-(n2_t2 - (i-j))+1);
            vxy_vec_Pwr = [fxy_vec_Pwr ; temp];
        end
    end
end

%% Get Coefficients of polynomial d
d_x_poly_pwr = BuildPoly_Pwr(roots_d_x);
d_y_poly_pwr = BuildPoly_Pwr(roots_d_y);
dxy_matrix_Pwr = d_x_poly_pwr * d_y_poly_pwr';
end

