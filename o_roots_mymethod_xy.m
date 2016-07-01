function [wx,wy,wxy] = o_roots_mymethod_xy(wx,wy,vDegt_wx,vDegt_wy)

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
        alpha = 1;
        th1 = 1;
        th2 = 1;
        
        [~,nCols] = size(fxy_matrix_n);
        [nRows,~] = size(gxy_matrix_n);
        t1 = nRows -1;
        t2 = nCols -1;
        
        
        % Get quotients
        [uxy_matrix_clc,vxy_matrix_clc] = GetQuotients(fxy_matrix_n,gxy_matrix_n,t1,t2,alpha,th1,th2);
        
        %[uxy_matrix_clc,vxy_matrix_clc] = GetQuotients_both(fxy_matrix_n,gxy_matrix_n,m,n,t,t1,t2,alpha,th1,th2);
        
        % Get the GCD dxy
        dxy_matrix_clc = GetGCDCoefficients(fxy_matrix_n,gxy_matrix_n,uxy_matrix_clc,vxy_matrix_clc,alpha, th1, th2);

%         m = vDegt_wx(i);
%         n = vDegt_wy(i);
%         lower_lim = 0;
%         upper_lim = min(m,n);
%         
%         [~,~,dxy_matrix_clc,~,~,t,t1,t2] = o_gcd_mymethod(fxy_matrix_n,gxy_matrix_n,m,n,[lower_lim, upper_lim]);
%         
        %Overwrite wx and wy with new values
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


