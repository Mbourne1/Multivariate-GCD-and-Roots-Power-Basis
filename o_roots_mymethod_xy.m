function [wx,wy,wxy] = o_roots_mymethod_xy(wx,wy,vDegt_x,vDegt_y)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Perform a series of GCD calculations on the w_{x,i}s

LineBreakLarge()

% Get number of w_{x,i}
[~,nEntries_wx] = size(wx);



for i = 1:1:nEntries_wx
    
    
    [~,cols] = size(wx{i});
    if cols > 1
        
        fxy_matrix_n = wx{i};
        gxy_matrix_n = wy{i};
        opt_alpha = 1;
        opt_theta_1 = 1;
        opt_theta_2 = 1;
        
        [~,t2] = GetDegree(fxy_matrix_n);
        [t1,~] = GetDegree(gxy_matrix_n);
        
        % Get quotients
        [uxy_matrix_clc,vxy_matrix_clc] = GetQuotients(fxy_matrix_n,gxy_matrix_n,t1,t2,opt_alpha,opt_theta_1,opt_theta_2);
        
        % Get the GCD dxy
        dxy_matrix_clc = GetGCDCoefficients(fxy_matrix_n,gxy_matrix_n,uxy_matrix_clc,vxy_matrix_clc,opt_alpha, opt_theta_1, opt_theta_2);
        
        % Overwrite wx and wy with new values
        wxy{i} = dxy_matrix_clc;
        wx{i} = uxy_matrix_clc;
        wy_new = Deconvolve_Bivariate(wy{i},dxy_matrix_clc);
        wy{i} = vxy_matrix_clc;
        
        fprintf([mfilename ' : ' sprintf('Roots of degree %i',i)])
        display(wxy{i})
        display(wx{i})
        display(wy{i})
        
    end
end


