function [] = o_Intersection_Implicit_Curves
% Given two implicit curves, calculate the points of intersection

global PLOT_GRAPHS 
PLOT_GRAPHS = 'n'

% C1 := f(x,y) = 0.
% C2 := g(x,y) = 0.

x = sym('x');
y = sym('y');

C1 = ...
    [
    -5 0 1;
    0 0 0;
    1 0 0;
    ];

C1_sym = x^2 + y^2 - 5 ;

C2 = ...
    [
        0   -1
        1   0
    ];

C2_sym = x-y;
C3_sym = C1_sym - C2_sym

figure('name','Plot')
hold on
ezplot(C1_sym)
ezplot(C2_sym)
ezplot(C3_sym)
hold off


C3 = ...
    [
    -5  1   1;
    -1  0   0;
    1   0   0;
    ];
m = 2;

o_roots_mymethod(C3,m)



end