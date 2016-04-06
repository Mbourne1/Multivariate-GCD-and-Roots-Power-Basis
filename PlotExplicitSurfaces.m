function [] = PlotExplicitSurfaces(aPoly)
% Given a set of Explicitly defined surfaces z=f(x,y), print the plots.

nPolys = length(aPoly);

% Given an array of bivariate polynomials of the form f(x,y), plot the
% surfaces z = f(x,y).

lwr_bound_x = -1;
upr_bound_x = 10;

lwr_bound_y = -1;
upr_bound_y = 10;

inc =0.1;

x_val_vec = lwr_bound_x : inc : upr_bound_x;
y_val_vec = lwr_bound_y : inc : upr_bound_y;

[p,q] = meshgrid(x_val_vec, y_val_vec);

figure('name',strcat(mfilename() , ' - Surface Plot'))
hold on
%zlim([-2,2])

for k = 1:1:nPolys
    
    % Initialise the matrix of z values
    z{k} = zeros(length(x_val_vec),length(y_val_vec));
    
    % Get the z values
    for i = 1:1:length(x_val_vec)
        for j = 1:1:length(y_val_vec)
            % Get the z values
            z{k}(i,j) = Evaluate_PowerPoly_Bivariate(x_val_vec(i),y_val_vec(j),...
                aPoly{k});
            
        end
    end
    
    % Produce surfaces
    s{k} = surf(p,q,z{k}');
    %set(s{k},'FaceColor',[1 0 0], 'FaceAlpha',0.5, 'EdgeAlpha', 0);
    % Set the colour of the surface using a binary representation of the
    % intger 'k'
    RGB_ = [bitget(k,1) bitget(k,2) bitget(k,3)];
    set(s{k},'FaceColor',RGB_, 'EdgeAlpha', 0);
    alpha(s{k},0.5);
end

legend(gca,'show')
hold off


end