function [] = o_IntersectionImplicitParametric()
% Get the intersection between an Implicitly defined surface and a
% parametrically defined surface.

% S_Explicit
% An Explicitly defined surface z = f(x,y)

% The surface S_{1} is defined by z = f(x,y)
S_Explicit = ...
    [
        0 1 
        1 0
    ];

s = sym('s');
t = sym('t');

S1_funx = s;
S1_funy = t;
S1_funz = s + t;




x = sym('x');
y = sym('y');
z = sym('z');

% Print the explicit surface
surfaces{1} = S_Explicit;



S1 = x + y - z;


%% 
% Get the parametrically defined surface
x_st = ...
    [
        5   0
        1   0
    ];
S2_funx = GetSymBivariatePoly(x_st)

y_st = ...
    [
        0   1
        0   0
    ];
S2_funy = GetSymBivariatePoly(y_st)

z_st = ...
    [
        0   1
        1   1
    ];
S2_funz = GetSymBivariatePoly(z_st)

%%
figure('name','Surface Plot')
hold on
Surf1 = ezsurf(S1_funx,S1_funy,S1_funz)
Surf2 = ezsurf(S2_funx,S2_funy,S2_funz)
hold off

f_st = subs(S1,{x y z},{x_st, y_st,z_st})

f_st = double(f_st);
[r,c] = size(f_st);

m1 = r - 1;
m2 = c - 1;
m = m1 + m2;



o_roots_mymethod(f_st,m)

% S_Parametric
% A Parametrically defined Surface S(x,y,z) = x(s,t), y(s,t), z(s,t)


end
