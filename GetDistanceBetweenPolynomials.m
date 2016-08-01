function [dist] = GetDistanceBetweenPolynomials(fxy_exact,fxy_computed,name)
% Given the two polynomials f(x) exact and f(x) computed. Calculate the
% distance between them.

LineBreakLarge()
display(fxy_exact)
display(fxy_computed)

% Normalise f(x,y) exact
fxy_exact = NormaliseMatrix(fxy_exact);

% Normalise f(x,y) computed
fxy_computed = NormaliseMatrix(fxy_computed);

DISTANCE_METRIC = 'Distance';


switch DISTANCE_METRIC
    case 'Angle'
        % Get angle between vectors
        vFxy_exact = GetAsVector(fxy_exact);
        vFxy_comp = GetAsVector(fxy_computed);
        
        dist = GetAngle(vFxy_exact,vFxy_comp);
        
        fprintf([mfilename ' : ' sprintf('Angle between %s exact and %s computed \t %2.4e \n',name,name,dist)])
        
    case 'Distance'
        
        dist = norm(fxy_exact - fxy_computed) ./ norm(fxy_exact);
        fprintf([mfilename ' : ' sprintf('Distance of %s exact from %s computed a - b / a : \t %2.4e \n',name,name,dist)]);
        
    otherwise
        error('err')
end



end

function costheta = GetAngle(a,b)
costheta = dot(a,b)/(norm(a)*norm(b));
end