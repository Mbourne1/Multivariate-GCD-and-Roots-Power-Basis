function [fxy_lr,gxy_lr,alpha_lr,theta1_lr,theta2_lr,x_lr] = ...
    SNTLN(fxy_n,gxy_n,alpha,th1,th2,m,n,t,t1,t2,opt_col)

% Initialise Global Settings
global SETTINGS

% Get method to compute SNTLN
switch SETTINGS.CALC_METHOD
    case 'Total'
        % 
        fprintf('No STLN method specified in terms of total degree \n')
        
    case 'Relative'
        % Get the SNTLN of the Sylvester matrix
        [fxy_lr,gxy_lr,alpha_lr,theta1_lr,theta2_lr,x_lr] = ...
            SNTLN_Respective(fxy_n,gxy_n,alpha,th1,th2,t1,t2,opt_col);
        
    case 'Both'
        % Get the SNTLN of the Sylvester matrix
        [fxy_lr,gxy_lr,alpha_lr,theta1_lr,theta2_lr,x_lr] = ...
            SNTLN_Both(fxy_n,gxy_n,alpha,th1,th2,m,n,t,t1,t2,opt_col);
end

end