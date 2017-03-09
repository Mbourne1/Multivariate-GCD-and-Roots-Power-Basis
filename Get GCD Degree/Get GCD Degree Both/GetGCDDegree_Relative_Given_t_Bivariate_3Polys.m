function [t1,t2] = GetGCDDegree_Relative_Given_t_Bivariate_3Polys(fxy, gxy, hxy, m, n, o, t, myLimits_t1, myLimits_t2, limits_t1, limits_t2)
% Get the degree structure (t_{1} and t_{2}) of the GCD d(x,y) of the two
% polynomials f(x,y) and g(x,y)
%
% % Inputs.
%
% [fxy, gxy, hxy] : Coefficient matrices of polynomials f(x,y), g(x,y) and
% h(x,y)
%
% [m, n, o] : Total degree of polynomial f(x,y), g(x,y) and h(x,y)
%
% t : Total degree of GCD d(x,y)
%
% % Outputs
%
% t1 : degree of d(x,y) with respect to x
% 
% t2 : degree of d(x,y) with respect to y

% Note this file differs from 'GetGCDDegreeRelative' since we make use of
% the computed value of 't' given t, for any k1,k2 pair, we know the
% location of certain zeros (if they exist in the lower right triangle) for
% polynomials u(x,y) and v(x,y). We remove the corresponding columns of
% S(f,g)

lowerLimit_t1 = myLimits_t1(1);
upperLimit_t1 = myLimits_t1(2);
lowerLimit_t2 = myLimits_t2(1);
upperLimit_t2 = myLimits_t2(2);

nSubresultants_k1 = upperLimit_t1 - lowerLimit_t1 + 1;
nSubresultants_k2 = upperLimit_t2 - lowerLimit_t2 + 1;

arr_SingularValues = cell( nSubresultants_k1, nSubresultants_k2);
arr_R = cell(nSubresultants_k1, nSubresultants_k2);
arr_R1 = cell( nSubresultants_k1, nSubresultants_k2);
arr_Diagonals_R1 = cell( nSubresultants_k1, nSubresultants_k2);


% For each of the pairs [k1,k2]
for i1 = 1:1: nSubresultants_k1
    for i2 = 1:1: nSubresultants_k2
        
        
        k1 = lowerLimit_t1 + (i1 - 1);
        k2 = lowerLimit_t2 + (i2 - 1);
        
        
        % Build the Sylvester matrix S_{k,k1,k2}
        Skk1k2 = BuildT_Both_Bivariate_3Polys(fxy, gxy, hxy, m, n, o, t, k1, k2);

        % Get the singular values of S_{k_{1},k_{2}}
        arr_SingularValues{i1,i2} = svd(Skk1k2);

    end
end



global SETTINGS
% R1 Row Norms
% R1 Row Diagonals
% Singular Values
% Residuals

switch SETTINGS.RANK_REVEALING_METRIC
    case 'Singular Values'
        mat_MinimumSingularValues = zeros( nSubresultants_k1, nSubresultants_k2);
        
        for i1 = 1:1: nSubresultants_k1
            for i2 = 1:1: nSubresultants_k2
                
                %k1 = lowerLimit_k1 + (i1 - 1);
                %k2 = lowerLimit_k2 + (i2 - 1);
                
                mat_MinimumSingularValues(i1,i2) = min(arr_SingularValues{i1,i2});
                
            end
        end
        
        mat_metric = mat_MinimumSingularValues;
        
        plotSingularValues_degreeRelative(arr_SingularValues, myLimits_t1, myLimits_t2, limits_t1, limits_t2);
        plotMinimumSingularValues_degreeRelative(mat_MinimumSingularValues, myLimits_t1, myLimits_t2, limits_t1, limits_t2);
        
        
    case 'R1 Row Norms'
        error('Not Completed')
    case 'R1 Row Diagonals'
        error('Not Completed')
    case 'Residuals'
        error('Not Completed')
    otherwise 
        error('Error : Not a valid rank revealing metric')
end


% Compute the degree of the GCD
delta_x = diff(log10(mat_metric),1,1);
vec_delta_x = sum(delta_x,2);
[~, idx] = max(vec_delta_x);
t1 = lowerLimit_t1 + idx - 1;

delta_y = diff(log10(mat_metric),1,2);
vec_delta_y = sum(delta_y,1);
[~, idx] = max(vec_delta_y);
t2 = lowerLimit_t2 + idx - 1;


fprintf([mfilename ' : ' 'The Calculated Degree of the GCD is given by \n'])
fprintf([mfilename ' : ' sprintf('Degree of GCD wrt x : t1 = %i\n',t1)])
fprintf([mfilename ' : ' sprintf('Degree of GCD wrt y : t2 = %i\n',t2)])
end
