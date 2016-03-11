function [] = 3DCurve_Fitting(m)

% Get a set of points in 3 dimensions

CP = ...
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
scatter3(X,Y,Z)
hold off

% Fit a polynomial curve


end



