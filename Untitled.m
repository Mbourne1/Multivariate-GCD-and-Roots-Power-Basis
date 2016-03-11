function [  ] = Untitled()
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% Get a set of points in 3 dimensions

CP = ...
    [
    1   2   3;
    1   7   9;
    1   9   12;
    ]

A = ...
    [
    1     1     1     1     1     1;
    4    14    49     2     7     1;
    9    27    81     3     9     1;
    ];

b = ...
    [
    1
    9
    12
    ];

% Plot the points
figure()
hold on
X = CP(1,:)
Y = CP(2,:)
Z = CP(3,:)
scatter3(X,Y,Z,'filled')
hold off

% Fit a polynomial curve
coef = pinv(transpose(A)*A)*transpose(A)*b

x = sym('x')
y = sym('y')
mypoly = -0.0616 - 0.2614*x + 0.1662*y +0.2787*x*x + 0.5544*x*y + 0.3237*y*y;

x_array = meshgrid(0:10,0:10)'
y_array = meshgrid(0:10,0:10)


for i = 0:1:10
    for j = 0:1:10
        z_array(i+1,j+1) = subs(mypoly,{x,y},{i,j});
    end
end

z_array

plot3(x_array,y_array,z_array)

end

